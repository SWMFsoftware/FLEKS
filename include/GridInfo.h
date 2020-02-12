#ifndef _GRIDINFO_H_
#define _GRIDINFO_H_

#include <cmath>
#include <iostream>

// Give it an unusual name so that it is less likely to be confilict with local
// variables.
extern int* _GRID_INFO_;

inline void reset_grid_info() {
  const int nx = _GRID_INFO_[1];
  const int ny = _GRID_INFO_[2];
  const int nz = _GRID_INFO_[3];
  const int nInt = 4 + ceil(nx * ny * nz / 8.0 / sizeof(int));
  for (int i = 4; i < nInt; i++) {
    _GRID_INFO_[i] = 0;
  }
}

inline void init_grid_info(int nPatchSize, int nx, int ny, int nz) {
  int nInt = 4 + ceil(nx * ny * nz / 8.0 / sizeof(int));
  _GRID_INFO_ = new int[nInt];
  _GRID_INFO_[0] = nPatchSize;
  _GRID_INFO_[1] = nx;
  _GRID_INFO_[2] = ny;
  _GRID_INFO_[3] = nz;
  reset_grid_info();
}

inline void set_point_status(int i, int j, int k, int val) {
  // const int nPatchSize = *(int*)_GRID_INFO_;
  const int nx = _GRID_INFO_[1];
  const int ny = _GRID_INFO_[2];
  const int nz = _GRID_INFO_[3];

  const int nBitInt = sizeof(int) * 8;
  const int nBitShift = nBitInt * 4 + i * (ny * nz) + j * (nz) + k;

  const int nIntShift = floor(nBitShift / nBitInt);

  const int nBitShiftLoc = nBitShift - nIntShift * nBitInt;

  // Note: i is a reference.
  int& number = _GRID_INFO_[nIntShift];

  if (val == 1) {
    number |= 1 << nBitShiftLoc;
  } else if (val == 0) {
    number &= ~(1 << nBitShiftLoc);
  } else {
    abort();
  }

  std::cout << "nx = " << nx << " ny = " << ny << " nz = " << nz
            << " sizeof(int) = " << sizeof(int) << " nBitShift = " << nBitShift
            << std::endl;
}

inline int get_point_status(int i, int j, int k) {
  const int nx = _GRID_INFO_[1];
  const int ny = _GRID_INFO_[2];
  const int nz = _GRID_INFO_[3];

  const int nBitInt = sizeof(int) * 8;
  const int nBitShift = nBitInt * 4 + i * (ny * nz) + j * (nz) + k;

  const int nIntShift = floor(nBitShift / nBitInt);

  const int nBitShiftLoc = nBitShift - nIntShift * nBitInt;

  int number = _GRID_INFO_[nIntShift];

  int val = (number >> nBitShiftLoc) & 1;
  return val;
}

inline void free_grid_info() { delete[] _GRID_INFO_; }

#endif