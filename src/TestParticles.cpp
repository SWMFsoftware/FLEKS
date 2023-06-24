#include <AMReX_Utility.H>

#include "TestParticles.h"

using namespace amrex;

TestParticles::TestParticles(amrex::AmrCore* amrcore,
                             FluidInterface* const fluidIn, TimeCtr* const tcIn,
                             const int speciesID, const amrex::Real charge,
                             const amrex::Real mass, int id)
    : Particles(amrcore, fluidIn, tcIn, speciesID, charge, mass,
                IntVect(AMREX_D_DECL(1, 1, 1))) {
  gridID = id;

  gridName = std::string("FLEKS") + std::to_string(gridID);
  printPrefix = gridName + ": ";

  for (int iDim = 0; iDim < nDim; iDim++)
    nPartPerCell[iDim] = 2;

  outputDir = component + "/plots/test_particles";

  nInitPart = 0;
}

//==========================================================
void TestParticles::move_and_save_particles(const amrex::MultiFab& nodeEMF,
                                            const amrex::MultiFab& nodeBMF,
                                            amrex::Real dt, amrex::Real dtNext,
                                            amrex::Real tNowSI, bool doSave) {
  if (is_neutral()) {
    move_and_save_neutrals(dt, tNowSI, doSave);
  } else {
    move_and_save_charged_particles(nodeEMF, nodeBMF, dt, dtNext, tNowSI,
                                    doSave);
  }
}

//==========================================================
void TestParticles::move_and_save_charged_particles(
    const amrex::MultiFab& nodeEMF, const amrex::MultiFab& nodeBMF,
    amrex::Real dt, amrex::Real dtNext, amrex::Real tNowSI, bool doSave) {
  timing_func("TestParticles::move_charged_particles");

  Real dtLoc = 0.5 * (dt + dtNext);

  if (particlePosition == NonStaggered) {
    update_position_to_half_stage(nodeEMF, nodeBMF, dt);
    dtLoc = 0.5 * dt;
  }

  const Real qdto2mc = charge / mass * 0.5 * dt;

  const int iLev = 0;
  for (ParticlesIter<nPTPartReal, nPTPartInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    const Array4<Real const>& nodeEArr = nodeEMF[pti].array();
    const Array4<Real const>& nodeBArr = nodeBMF[pti].array();

    const Array4<int const>& status = cellStatus[iLev][pti].array();
    // cellStatus[iLev][pti] is a FAB, and the box returned from the box()
    // method already contains the ghost cells.
    const Box& bx = cellStatus[iLev][pti].box();
    const IntVect lowCorner = bx.smallEnd();
    const IntVect highCorner = bx.bigEnd();

    auto& particles = pti.GetArrayOfStructs();
    for (auto& p : particles) {
      if (p.idata(iRecordCount_) >= nPTRecord) {
        Abort("Error: there is not enough allocated memory to store the "
              "particle record!!");
      }

      Real up = p.rdata(iup_);
      Real vp = p.rdata(ivp_);
      Real wp = p.rdata(iwp_);
      const Real xp = p.pos(ix_);
      const Real yp = p.pos(iy_);
      const Real zp = p.pos(iz_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      RealVect dShift;
      for (int i = 0; i < 3; i++) {
        dShift[i] = (p.pos(i) - plo[iLev][i]) * invDx[iLev][i];
        loIdx[i] = fastfloor(dShift[i]);
        dShift[i] = dShift[i] - loIdx[i];
      }

      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      Real Bxl = 0, Byl = 0, Bzl = 0; // should be bp[3];
      Real Exl = 0, Eyl = 0, Ezl = 0;
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++) {
            const int iNodeX = loIdx[ix_] + ii;
            const int iNodeY = loIdx[iy_] + jj;
            const int iNodeZ = loIdx[iz_] + kk;
            const Real& c0 = coef[ii][jj][kk];
            Bxl += nodeBArr(iNodeX, iNodeY, iNodeZ, ix_) * c0;
            Byl += nodeBArr(iNodeX, iNodeY, iNodeZ, iy_) * c0;
            Bzl += nodeBArr(iNodeX, iNodeY, iNodeZ, iz_) * c0;

            Exl += nodeEArr(iNodeX, iNodeY, iNodeZ, ix_) * c0;
            Eyl += nodeEArr(iNodeX, iNodeY, iNodeZ, iy_) * c0;
            Ezl += nodeEArr(iNodeX, iNodeY, iNodeZ, iz_) * c0;
          }

      Real gamma = 1;
      Real invGamma = 1. / gamma;

      if (isRelativistic) {
        // Convert: vel -> gamma*vel
        const Real v2 = up * up + vp * vp + wp * wp;
        if (v2 > 1) {
          Abort("Error: particle speed is fast than light!");
        }
        invGamma = sqrt(1 - v2);
        gamma = 1 / invGamma;
        up *= gamma;
        vp *= gamma;
        wp *= gamma;
      }

      // Half step acceleration
      const Real ut = up + qdto2mc * Exl;
      const Real vt = vp + qdto2mc * Eyl;
      const Real wt = wp + qdto2mc * Ezl;

      if (isRelativistic) {
        const Real p2 = ut * ut + vt * vt + wt * wt;
        gamma = sqrt(1 + p2);
        invGamma = 1. / gamma;
      }

      const double Omx = qdto2mc * Bxl * invGamma;
      const double Omy = qdto2mc * Byl * invGamma;
      const double Omz = qdto2mc * Bzl * invGamma;

      // end interpolation
      const Real omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const Real denom = 1.0 / (1.0 + omsq);

      // const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
      const Real udotOm = ut * Omx + vt * Omy + wt * Omz;
      // solve the velocity equation
      const Real uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
      const Real vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
      const Real wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;

      double unp1 = 2.0 * uavg - up;
      double vnp1 = 2.0 * vavg - vp;
      double wnp1 = 2.0 * wavg - wp;

      if (isRelativistic) {
        // Convert: gamma*vel -> vel
        const Real p2 = unp1 * unp1 + vnp1 * vnp1 + wnp1 * wnp1;
        gamma = sqrt(1 + p2);
        invGamma = 1. / gamma;

        unp1 *= invGamma;
        vnp1 *= invGamma;
        wnp1 *= invGamma;
      }

      p.rdata(ix_) = unp1;
      p.rdata(iy_) = vnp1;
      p.rdata(iz_) = wnp1;

      p.pos(ix_) = xp + unp1 * dtLoc;
      p.pos(iy_) = yp + vnp1 * dtLoc;
      p.pos(iz_) = zp + wnp1 * dtLoc;

      if (doSave) {
        const int i0 = record_var_index(p.idata(iRecordCount_));
        p.rdata(i0 + iTPt_) = tNowSI;
        p.rdata(i0 + iTPu_) = unp1;
        p.rdata(i0 + iTPv_) = vnp1;
        p.rdata(i0 + iTPw_) = wnp1;
        p.rdata(i0 + iTPx_) = xp + unp1 * 0.5 * dt;
        p.rdata(i0 + iTPy_) = yp + vnp1 * 0.5 * dt;
        p.rdata(i0 + iTPz_) = zp + wnp1 * 0.5 * dt;

        if (ptRecordSize > iTPBx_) {
          p.rdata(i0 + iTPBx_) = Bxl;
          p.rdata(i0 + iTPBy_) = Byl;
          p.rdata(i0 + iTPBz_) = Bzl;
        }

        if (ptRecordSize > iTPEx_) {
          p.rdata(i0 + iTPEx_) = Exl;
          p.rdata(i0 + iTPEy_) = Eyl;
          p.rdata(i0 + iTPEz_) = Ezl;
        }

        p.idata(iRecordCount_) = p.idata(iRecordCount_) + 1;
      }
      // Print() << "p = " << p << std::endl;

      // Mark for deletion
      if (is_outside_ba(p, status, lowCorner, highCorner)) {
        p.id() = -1;
      }

    } // for p
  }   // for pti

  // This function distributes particles to proper processors and apply
  // periodic boundary conditions if needed.
  Redistribute();
}

//==========================================================
void TestParticles::move_and_save_neutrals(amrex::Real dt, amrex::Real tNowSI,
                                           bool doSave) {
  timing_func("TestParticles::move_neutrals");

  const int iLev = 0;
  for (ParticlesIter<nPTPartReal, nPTPartInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    const Array4<int const>& status = cellStatus[iLev][pti].array();
    // cellStatus[iLev][pti] is a FAB, and the box returned from the box()
    // method already contains the ghost cells.
    const Box& bx = cellStatus[iLev][pti].box();
    const IntVect lowCorner = bx.smallEnd();
    const IntVect highCorner = bx.bigEnd();

    auto& particles = pti.GetArrayOfStructs();
    for (auto& p : particles) {
      if (p.idata(iRecordCount_) >= nPTRecord) {
        Abort("Error: there is not enough allocated memory to store the "
              "particle record!!");
      }

      Real up = p.rdata(iup_);
      Real vp = p.rdata(ivp_);
      Real wp = p.rdata(iwp_);
      const Real xp = p.pos(ix_);
      const Real yp = p.pos(iy_);
      const Real zp = p.pos(iz_);

      p.pos(ix_) = xp + up * dt;
      p.pos(iy_) = yp + vp * dt;
      p.pos(iz_) = zp + wp * dt;

      if (doSave) {
        const int i0 = record_var_index(p.idata(iRecordCount_));
        p.rdata(i0 + iTPt_) = tNowSI;
        p.rdata(i0 + iTPu_) = p.rdata(iup_);
        p.rdata(i0 + iTPv_) = p.rdata(ivp_);
        p.rdata(i0 + iTPw_) = p.rdata(iwp_);
        p.rdata(i0 + iTPx_) = p.pos(ix_);
        p.rdata(i0 + iTPy_) = p.pos(iy_);
        p.rdata(i0 + iTPz_) = p.pos(iz_);

        p.idata(iRecordCount_) = p.idata(iRecordCount_) + 1;
      }

      // Mark for deletion
      if (is_outside_ba(p, status, lowCorner, highCorner)) {
        p.id() = -1;
      }
    }
  }

  // This function distributes particles to proper processors and apply
  // periodic boundary conditions if needed.
  Redistribute();
}

//======================================================================
void TestParticles::read_test_particle_list(
    const amrex::Vector<std::string>& listFiles) {
  std::string funcName = "TP::read_test_particle_list";
  timing_func(funcName);

  std::string listName;
  for (int i = 0; i < listFiles.size(); i++) {
    // Assume the particle list file's name is *N.dat, where N is
    // the species ID (Usually 0 or 1).
    if (listFiles[i].find(std::to_string(speciesID) + ".dat") !=
        std::string::npos) {
      listName = listFiles[i];
      break;
    }
  }

  if (listName.empty())
    return;

  vIDs.clear();

  // Read particle IDs
  std::ifstream source;
  source.open(listName, std::ifstream::in);
  PID id;
  while (source >> id.cpu >> id.id) {
    vIDs.push_back(id);
    id.flag = true;
  }
  source.close();

  std::sort(vIDs.begin(), vIDs.end());
}

//======================================================================
// Trace the PIC particles that are in the list (listFiles).
void TestParticles::add_test_particles_from_pic(Particles<>* pts) {
  std::string funcName = "TP::add_test_particles_from_pic";
  timing_func(funcName);

  if (vIDs.size() == 0)
    return;

  const int iLev = 0;

  const auto& lOther = pts->GetParticles(iLev);

  // Loop through the PIC particles.
  for (MFIter mfi = pts->MakeMFIter(iLev); mfi.isValid(); ++mfi) {
    auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());

    if (lOther.find(index) == lOther.end())
      continue;

    const auto& tileOther = lOther.at(index);

    if (tileOther.numParticles() == 0)
      continue;

    auto& particles = GetParticles(iLev)[index];

    const auto& aosOther = tileOther.GetArrayOfStructs();

    for (auto pOther : aosOther) {
      PID id;
      id.cpu = pOther.cpu();
      id.id = pOther.id();

      auto it = std::lower_bound(vIDs.begin(), vIDs.end(), id);

      // This particle is on the test particle list. Generate a new test
      // particle.
      if (*it == id) {
        ParticleType p;
        p.cpu() = pOther.cpu();
        p.id() = pOther.id();
        p.pos(ix_) = pOther.pos(ix_);
        p.pos(iy_) = pOther.pos(iy_);
        p.pos(iz_) = pOther.pos(iz_);
        p.rdata(iup_) = pOther.rdata(iup_);
        p.rdata(ivp_) = pOther.rdata(ivp_);
        p.rdata(iwp_) = pOther.rdata(iwp_);
        p.rdata(iqp_) = pOther.rdata(iqp_);

        p.idata(iRecordCount_) = 0;

        particles.push_back(p);

        // Remove this id from the test particle list later.
        (*it).flag = false;
      }
    }
  }

  { // If a test particle has been generated on one proc, all proc should remove
    // this particle from the test particle list (vIDs). Use a bit array to do
    // mpi_allreduce to improve performance.
    BitArray ba(vIDs.size());
    for (unsigned int i = 0; i < vIDs.size(); i++) {
      if (vIDs[i].flag) {
        ba.set(i, 1);
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, ba.get(), ba.size_int(), MPI_INT, MPI_BAND,
                  amrex::ParallelDescriptor::Communicator());

    for (unsigned int i = 0; i < vIDs.size(); i++) {
      if (ba.get(i) == 0) {
        vIDs[i].flag = false;
      }
    }
  }

  // After move, vIDs is empty.
  std::vector<PID> tmp = std::move(vIDs);

  for (auto t : tmp) {
    if (t.flag)
      vIDs.push_back(t);
  }
}

//======================================================================
void TestParticles::add_test_particles_from_fluid(Vector<Vel> tpStates) {
  std::string funcName = "TP::add_test_particles_from_fluid";
  timing_func(funcName);
  Print() << funcName << " : nInitPart = " << nInitPart
          << " : current number = " << TotalNumberOfParticles(true, false)
          << std::endl;

  const int iLev = 0;

  Vel tpVel;
  for (auto& state : tpStates) {
    if (state.tag == speciesID) {
      tpVel = state;
    }
  }

  for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
    const auto& status = cellStatus[iLev][mfi].array();
    const Box& bx = mfi.validbox();
    const IntVect lo = IntVect(bx.loVect());
    const IntVect hi = IntVect(bx.hiVect());

    IntVect idxMin = lo, idxMax = hi;

    for (int i = idxMin[ix_]; i <= idxMax[ix_]; i += nIntervalCell[ix_])
      for (int j = idxMin[iy_]; j <= idxMax[iy_]; j += nIntervalCell[iy_])
        for (int k = idxMin[iz_]; k <= idxMax[iz_]; k += nIntervalCell[iz_])
          if (iPartRegion == iRegionUniform_ ||
              (iPartRegion == iRegionBoundary_ &&
               bit::is_edge(status(i, j, k)))) {
            add_particles_cell(iLev, mfi, i, j, k, *fi,
                               { AMREX_D_DECL(0, 0, 0) }, tpVel);
          }
  }
}

//======================================================================
bool TestParticles::write_particles(int cycle) {
  std::string funcName = "TP::write_particles";
  timing_func(funcName);

  int nPartLoc = TotalNumberOfParticles(false, true);
  int nPartAhead = 0;

  gather_accumulate_and_scatter(nPartLoc, nPartAhead);

  unsigned long long int nByteLoc = loop_particles("count_record_size");
  unsigned long long int nByteAhead = 0;
  gather_accumulate_and_scatter(nByteLoc, nByteAhead);

  Vector<char> dataBuffer;
  dataBuffer.resize(nByteLoc);
  loop_particles("copy_record", dataBuffer.data(), dataBuffer.size());

  // cpu + id + loc;
  const int listUnitSize = 2 * sizeof(int) + sizeof(unsigned long long);
  Vector<char> partList;
  partList.resize(nPartLoc * listUnitSize);
  loop_particles("get_record_loc", partList.data(), partList.size(),
                 nByteAhead);

  std::stringstream ss;
  ss << "n" << std::setfill('0') << std::setw(8) << cycle;

  amrex::UtilCreateDirectory(outputDir, 0755);
  std::string fileNamePartList = outputDir + "/" + gridName +
                                 "_particle_list_species_" +
                                 std::to_string(speciesID) + "_" + ss.str();

  std::string fileNamePartRecord = outputDir + "/" + gridName +
                                   "_particle_species_" +
                                   std::to_string(speciesID) + "_" + ss.str();

  Print() << printPrefix << "Saving test particles into " << fileNamePartRecord
          << std::endl;

  MPI_File recordFile, listFile;
  MPI_Status status;

  MPI_File_open(ParallelDescriptor::Communicator(), fileNamePartList.c_str(),
                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &listFile);
  MPI_File_write_at(listFile, nPartAhead * listUnitSize, partList.data(),
                    partList.size(), MPI_CHAR, &status);
  MPI_File_close(&listFile);

  MPI_File_open(ParallelDescriptor::Communicator(), fileNamePartRecord.c_str(),
                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &recordFile);
  MPI_File_write_at(recordFile, nByteAhead, dataBuffer.data(),
                    dataBuffer.size(), MPI_CHAR, &status);
  MPI_File_close(&recordFile);

  loop_particles("reset_record_counter", partList.data(), partList.size(),
                 nByteAhead);

  {
    if (ParallelDescriptor::IOProcessor()) {
      std::string headerName = outputDir + "/Header";

      std::ofstream headerFile;
      headerFile.open(headerName.c_str(), std::ios::out | std::ios::trunc);

      if (!headerFile.good())
        amrex::FileOpenFailed(headerName);

      headerFile << ptRecordSize << "\n";
    }
  }

  return true;
}

unsigned long long int TestParticles::loop_particles(
    std::string action, char* buff, unsigned long long int sizeLimit,
    unsigned long long int shift) {

  unsigned long long int nByteCount = 0;

  bool doCopyData = (action == "copy_record");
  bool doGetLoc = (action == "get_record_loc");
  bool doResetRecordCounter = (action == "reset_record_counter");

  int iPartCount = 0;
  const int listUnitSize = 2 * sizeof(int) + sizeof(unsigned long long);

  const int iLev = 0;
  for (ParticlesIter<nPTPartReal, nPTPartInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    auto& particles = pti.GetArrayOfStructs();
    for (auto& p : particles) {
      int nRecord = p.idata(iRecordCount_);

      // int: cpu + id + nRecord
      // Real: weight + ptRecordSize*nRecord
      int nBytePerPart =
          3 * sizeof(int) + (1 + ptRecordSize * nRecord) * sizeof(float);

      if (doResetRecordCounter) {
        p.idata(iRecordCount_) = 0;
      }

      if (doGetLoc) {
        int i2[2];
        i2[0] = p.cpu();
        i2[1] = p.id();
        memcpy(buff + iPartCount * listUnitSize, i2, 2 * sizeof(int));

        unsigned long long tmp = nByteCount + shift;
        memcpy(buff + iPartCount * listUnitSize + 2 * sizeof(int), &tmp,
               sizeof(unsigned long long));
      }

      if (doCopyData) {
        int iCountLoc = 0, sizeLoc = 0;
        int i3[3];
        i3[0] = p.cpu();
        i3[1] = p.id();
        i3[2] = nRecord;

        // cpu + id + nRecord
        sizeLoc = 3 * sizeof(int);
        memcpy(buff + nByteCount + iCountLoc, i3, sizeLoc);
        iCountLoc += sizeLoc;

        // weight
        sizeLoc = sizeof(float);
        float weight = (float)(p.rdata(iqp_) / charge * mass * no2outM);
        memcpy(buff + nByteCount + iCountLoc, &weight, sizeLoc);
        iCountLoc += sizeLoc;

        sizeLoc = sizeof(float) * ptRecordSize;
        for (int i = 0; i < nRecord; i++) {
          float recordData[ptRecordSize];

          const int i0 = record_var_index(i);

          // time is already in SI unit
          recordData[iTPt_] = (float)p.rdata(i0 + iTPt_);
          recordData[iTPx_] = (float)(p.rdata(i0 + iTPx_) * no2outL);
          recordData[iTPy_] = (float)(p.rdata(i0 + iTPy_) * no2outL);
          recordData[iTPz_] = (float)(p.rdata(i0 + iTPz_) * no2outL);
          recordData[iTPu_] = (float)(p.rdata(i0 + iTPu_) * no2outV);
          recordData[iTPv_] = (float)(p.rdata(i0 + iTPv_) * no2outV);
          recordData[iTPw_] = (float)(p.rdata(i0 + iTPw_) * no2outV);

          if (ptRecordSize > iTPBx_) {
            recordData[iTPBx_] = (float)(p.rdata(i0 + iTPBx_) * no2outB);
            recordData[iTPBy_] = (float)(p.rdata(i0 + iTPBy_) * no2outB);
            recordData[iTPBz_] = (float)(p.rdata(i0 + iTPBz_) * no2outB);
          }

          if (ptRecordSize > iTPEx_) {
            recordData[iTPEx_] = (float)(p.rdata(i0 + iTPEx_) * no2outE);
            recordData[iTPEy_] = (float)(p.rdata(i0 + iTPEy_) * no2outE);
            recordData[iTPEz_] = (float)(p.rdata(i0 + iTPEz_) * no2outE);
          }

          memcpy(buff + nByteCount + iCountLoc, &recordData, sizeLoc);
          iCountLoc += sizeLoc;
        }

        if (nByteCount + iCountLoc > sizeLimit) {
          AllPrint() << "nByteCount = " << nByteCount
                     << " iCountLoc = " << iCountLoc
                     << " sizeLimit = " << sizeLimit << std::endl;
          Abort("Error: memory may leaked!!");
        }
      }

      iPartCount++;
      nByteCount += nBytePerPart;
    }
  }
  return nByteCount;
}

void TestParticles::print_record_buffer(char* buffer,
                                        unsigned long long int nBuffer) {
  unsigned long long int count = 0;

  while (count < nBuffer) {
    AllPrint() << "-----Record-----------" << std::endl;
    AllPrint() << "Record: cpu = " << (*(int*)(buffer + count));
    count += sizeof(int);
    AllPrint() << " id = " << (*(int*)(buffer + count));
    count += sizeof(int);

    int nRecord = (*(int*)(buffer + count));
    AllPrint() << " nRecord = " << nRecord;
    count += sizeof(int);

    AllPrint() << " weight = " << (*(float*)(buffer + count)) << std::endl;
    count += sizeof(float);

    for (int i = 0; i < nRecord; i++) {
      for (int ii = 0; ii < ptRecordSize; ii++) {
        AllPrint() << " data" << ii << " = " << (*(float*)(buffer + count));
        count += sizeof(float);
      }
      AllPrint() << std::endl;
    }
    AllPrint() << "----------------" << std::endl;
  }
}

void TestParticles::reset_record_counter() {
  const int iLev = 0;
  for (ParticlesIter<nPTPartReal, nPTPartInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    auto& particles = pti.GetArrayOfStructs();
    for (auto& p : particles) {
      p.idata(iRecordCount_) = 0;
    }
  }
}
