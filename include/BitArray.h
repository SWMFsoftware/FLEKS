#ifndef _BITARRAY_H_
#define _BITARRAY_H_

#include <cmath>
#include <stdio.h>

class BitArray {

private:
  int* data;
  int nInt, nBit;

public:
  BitArray() {
    data = nullptr;
    nInt = 0;
    nBit = 0;
  }

  BitArray(int nBitIn) { init(nBitIn); }

  ~BitArray() { clear(); }

  int size_int() { return nInt; }
  int size_bit() { return nBit; }

  //====================================================
  void clear() {
    if (data != nullptr) {
      delete[] data;
    }
    nInt = 0;
    nBit = 0;
  }

  //====================================================
  void init(int nBitIn) {
    nBit = nBitIn;
    nInt = ceil(nBit / 8.0 / sizeof(int));

    data = new int[nInt];

    for (int i = 0; i < nInt; ++i) {
      data[i] = 0;
    }
  }

  //====================================================
  void set(int i, int val) {
    if (i > nBit)
      abort();

    const int intSize = sizeof(int) * 8;

    const int iInt = floor(i / intSize);

    const int iBit = i % intSize;

    int& number = data[iInt];

    if (val == 1) {
      number |= 1 << iBit;
    } else if (val == 0) {
      number &= ~(1 << iBit);
    } else {
      abort();
    }
  }

  //====================================================
  int get(int i) {
    if (i > nBit)
      abort();

    const int intSize = sizeof(int) * 8;

    const int iInt = floor(i / intSize);

    const int iBit = i % intSize;

    int number = data[iInt];

    int val = (number >> iBit) & 1;

    return val;
  }

  //====================================================
  // Return an interger pointer
  int* get() { return data; }

  void print() {
    for (int i = 0; i < nBit; ++i) {
      printf("i = %i, bit = %i \n", i, get(i));
    }
  }
};

#endif