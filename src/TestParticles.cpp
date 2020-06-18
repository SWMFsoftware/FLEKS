#include <AMReX_Utility.H>

#include "TestParticles.h"

using namespace amrex;

TestParticles::TestParticles(const amrex::BoxArray& regionBAIn,
                             const amrex::Geometry& geom,
                             const amrex::DistributionMapping& dm,
                             const amrex::BoxArray& ba,
                             FluidInterface* const fluidIn, TimeCtr* const tcIn,
                             const int speciesID, const amrex::Real charge,
                             const amrex::Real mass, int domainIDIn)
    : Particles(regionBAIn, geom, dm, ba, fluidIn, tcIn, speciesID, charge,
                mass, IntVect(1, 1, 1)) {
  domainID = domainIDIn;

  {
    std::stringstream ss;
    ss << "FLEKS" << domainID;
    domainName = ss.str();
    printPrefix = domainName + ": ";
  }

  iFileCount = 0;
  iStep = 0;
  outputDir = "PC/plots/test_particles";

  nInitPart = 0;
}

void TestParticles::move_and_save_particles(const amrex::MultiFab& nodeEMF,
                                            const amrex::MultiFab& nodeBMF,
                                            amrex::Real dt, amrex::Real dtNext,
                                            amrex::Real tNowSI) {
  timing_func("TestParticles::mover");

  Real dtLoc = 0.5 * (dt + dtNext);
  const auto& plo = Geom(0).ProbLo();

  const Real qdto2mc = charge / mass * 0.5 * dt;

  const int lev = 0;
  for (ParticlesIter<nPTPartReal, nPTPartInt> pti(*this, lev); pti.isValid();
       ++pti) {
    const Array4<Real const>& nodeEArr = nodeEMF[pti].array();
    const Array4<Real const>& nodeBArr = nodeBMF[pti].array();

    const Array4<int const>& status = cellStatus[pti].array();
    // cellStatus[pti] is a FAB, and the box returned from the box() method
    // already contains the ghost cells.
    const Box& bx = cellStatus[pti].box();
    const IntVect lowCorner = bx.smallEnd();
    const IntVect highCorner = bx.bigEnd();

    auto& particles = pti.GetArrayOfStructs();
    for (auto& p : particles) {
      if (p.idata(iRecordCount_) >= nPTRecord) {
        Abort("Error: there is not enough allocated memory to store the "
              "particle record!!");
      }

      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);
      const Real xp = p.pos(ix_);
      const Real yp = p.pos(iy_);
      const Real zp = p.pos(iz_);

      //-----calculate interpolate coef begin-------------
      int loIdx[3];
      Real dShift[3];
      for (int i = 0; i < 3; i++) {
        dShift[i] = (p.pos(i) - plo[i]) * invDx[i];
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

      const double Omx = qdto2mc * Bxl;
      const double Omy = qdto2mc * Byl;
      const double Omz = qdto2mc * Bzl;

      // end interpolation
      const Real omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const Real denom = 1.0 / (1.0 + omsq);
      // solve the position equation
      const Real ut = up + qdto2mc * Exl;
      const Real vt = vp + qdto2mc * Eyl;
      const Real wt = wp + qdto2mc * Ezl;
      // const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
      const Real udotOm = ut * Omx + vt * Omy + wt * Omz;
      // solve the velocity equation
      const Real uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
      const Real vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
      const Real wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;

      const double unp1 = 2.0 * uavg - up;
      const double vnp1 = 2.0 * vavg - vp;
      const double wnp1 = 2.0 * wavg - wp;

      p.rdata(ix_) = unp1;
      p.rdata(iy_) = vnp1;
      p.rdata(iz_) = wnp1;

      p.pos(ix_) = xp + unp1 * dtLoc;
      p.pos(iy_) = yp + vnp1 * dtLoc;
      p.pos(iz_) = zp + wnp1 * dtLoc;

      const int iStart_ = record_var_index(p.idata(iRecordCount_));
      p.rdata(iStart_ + iRecordt_) = tNowSI;
      p.rdata(iStart_ + iRecordu_) = unp1;
      p.rdata(iStart_ + iRecordv_) = vnp1;
      p.rdata(iStart_ + iRecordw_) = wnp1;
      p.rdata(iStart_ + iRecordx_) = xp + unp1 * 0.5 * dt;
      p.rdata(iStart_ + iRecordy_) = yp + vnp1 * 0.5 * dt;
      p.rdata(iStart_ + iRecordz_) = zp + wnp1 * 0.5 * dt;

      p.idata(iRecordCount_) = p.idata(iRecordCount_) + 1;

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

  iStep++;
}

void TestParticles::add_test_particles(const iMultiFab& cellStatus) {
  std::string funcName = "TestParticles::add_test_particles";
  timing_func(funcName);
  Print() << funcName << " : nInitPart = " << nInitPart
          << " : current number = " << TotalNumberOfParticles(true, false) << std::endl;

  const int lev = 0;

  const Real partNumCtr = 0.2;

  for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const auto& status = cellStatus[mfi].array();
    const Box& bx = mfi.validbox();
    const IntVect lo = IntVect(bx.loVect());
    const IntVect hi = IntVect(bx.hiVect());

    IntVect idxMin = lo, idxMax = hi;

    for (int i = idxMin[ix_]; i <= idxMax[ix_]; ++i)
      for (int j = idxMin[iy_]; j <= idxMax[iy_]; ++j)
        for (int k = idxMin[iz_]; k <= idxMax[iz_]; ++k) {
          if (status(i, j, k) == iAddPTParticle_ && randNum() < partNumCtr) {
            add_particles_cell(mfi, i, j, k);
          }
        }
  }
}

//======================================================================
bool TestParticles::write_particles(bool forceOutput) {

  if (iStep % nPTRecord != 0 && !forceOutput)
    return false;

  const int nProc = ParallelDescriptor::NProcs();
  int nPartLoc = TotalNumberOfParticles(false, true);
  int nPartAhead = 0;

  gather_accumulate_and_scatter(nPartLoc, nPartAhead);

  unsigned long long int nByteLoc = loop_particles("count_record_size");
  unsigned long long int nByteAhead = 0;
  gather_accumulate_and_scatter(nByteLoc, nByteAhead);

  Vector<char> dataBuffer;
  dataBuffer.resize(nByteLoc);
  unsigned long long tmp;
  tmp = loop_particles("copy_record", dataBuffer.data(), dataBuffer.size());

  // cpu + id + loc;
  const int listUnitSize = 2 * sizeof(int) + sizeof(unsigned long long);
  Vector<char> partList;
  partList.resize(nPartLoc * listUnitSize);
  tmp = loop_particles("get_record_loc", partList.data(), partList.size(),
                       nByteAhead);

  std::stringstream ss;
  ss << std::setfill('0') << std::setw(4) << std::to_string(iFileCount);

  amrex::UtilCreateDirectory(outputDir, 0755);
  std::string fileNamePartList = outputDir + "/" + domainName +
                                 "_particle_list_species_" +
                                 std::to_string(speciesID) + "_" + ss.str();

  std::string fileNamePartRecord = outputDir + "/" + domainName +
                                   "_particle_species_" +
                                   std::to_string(speciesID) + "_" + ss.str();
  iFileCount++;

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

  return true;
}

unsigned long long int TestParticles::loop_particles(
    std::string action, char* buff, int sizeLimit,
    unsigned long long int shift) {

  unsigned long long int nByteCount = 0;

  bool doCopyData = (action == "copy_record");
  bool doGetLoc = (action == "get_record_loc");
  bool doResetRecordCounter = (action == "reset_record_counter");

  int iPartCount = 0;
  const int listUnitSize = 2 * sizeof(int) + sizeof(unsigned long long);

  const int lev = 0;
  for (ParticlesIter<nPTPartReal, nPTPartInt> pti(*this, lev); pti.isValid();
       ++pti) {
    auto& particles = pti.GetArrayOfStructs();
    for (auto& p : particles) {
      int nRecord = p.idata(iRecordCount_);

      // int: cpu + id + nRecord
      // Real: weight + (t, x, y, z, ux, uy, uz)*nRecord
      int nBytePerPart = 3 * sizeof(int) + (1 + 7 * nRecord) * sizeof(float);

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

        // t, x, y, z, u, v, w
        sizeLoc = sizeof(float) * 7;
        for (int i = 0; i < nRecord; i++) {
          float recordData[7];

          const int iStart_ = record_var_index(i);

          // time is already in SI unit
          recordData[iRecordt_] = (float)p.rdata(iStart_ + iRecordt_);
          recordData[iRecordx_] =
              (float)(p.rdata(iStart_ + iRecordx_) * no2outL);
          recordData[iRecordy_] =
              (float)(p.rdata(iStart_ + iRecordy_) * no2outL);
          recordData[iRecordz_] =
              (float)(p.rdata(iStart_ + iRecordz_) * no2outL);
          recordData[iRecordu_] =
              (float)(p.rdata(iStart_ + iRecordu_) * no2outV);
          recordData[iRecordv_] =
              (float)(p.rdata(iStart_ + iRecordv_) * no2outV);
          recordData[iRecordw_] =
              (float)(p.rdata(iStart_ + iRecordw_) * no2outV);

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
  int count = 0;

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
      for (int ii = 0; ii < 7; ii++) {
        AllPrint() << " data" << ii << " = " << (*(float*)(buffer + count));
        count += sizeof(float);
      }
      AllPrint() << std::endl;
    }
    AllPrint() << "----------------" << std::endl;
  }
}
