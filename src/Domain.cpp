#include "Domain.h"

void Domain::init() {
  {
    //---- Geometry initialization -------
    nGst = 2;

    nCell[ix_] = 32;
    nCell[iy_] = 32;
    nCell[iz_] = 32;

    nCellBlockMax = 16;

    for (auto& x : periodicity)
      x = 1;

    for (int i = 0; i < nDim; i++) {
      domainLo[i] = 0;
      domainHi[i] = nCell[i] - 1;

      domainRange.setLo(i, -1);
      domainRange.setHi(i, 1);
    }

    domainBox.setSmall(domainLo);
    domainBox.setBig(domainHi);

    coord = 0; // Cartesian

    geom.define(domainBox, &domainRange, coord, periodicity);

    ba.define(domainBox);
    ba.maxSize(nCellBlockMax);

    dm.define(ba);
  }

  {
    // EM field
    E.define(ba, dm, 3, nGst);
    E.setVal(0.0);
    B.define(ba, dm, 3, nGst);
    B.setVal(0.0);
  }

  {
    // Plasma
    nSpecies = 2;
    plasma.resize(nSpecies);
    for (auto& pl : plasma) {
      pl.define(ba, dm, nMoments, nGst);
      pl.setVal(0.0);
    }
  }

  AllPrint() << "Domain::init finished! "
             << "\n";
}
