#include "simdata.hpp"

namespace GFlowSimulation {

  SimData::SimData(GFlow *gflow) : Base(gflow), number(0), numberG(0), x(nullptr), v(nullptr), f(nullptr), 
                                  sg(nullptr), im(nullptr), type(nullptr), dataF(nullptr), dataI(nullptr),
                                  size(0) {};

  SimData::~SimData() {
    clean();
  }

  void SimData::clean() {
    int i;
    // Theses are arrays of arrays
    if (x) {
      for (i=0; i<DIMENSIONS; ++i) 
        if(x[i]) delete [] x[i];
      delete [] x;
      x = nullptr;
    }
    if (v) {
      for (i=0; i<DIMENSIONS; ++i) 
        if(v[i]) delete [] v[i];
      delete [] v;
      v = nullptr;
    }
    if (f) {
      for (i=0; i<DIMENSIONS; ++i) 
        if(f[i]) delete [] f[i];
      delete [] f;
      f = nullptr;
    }
    // These are just arrays
    if (sg) {
      delete [] sg;
      sg = nullptr;
    }
    if (im) {
      delete [] im;
      im = nullptr;
    }
    if (type) {
      delete [] type;
      type = nullptr;
    }
    // Delete additional data - these are arrays of arrays
    if (dataF) {
      for (i=0; i<DIMENSIONS; ++i) 
        if(dataF[i]) delete [] dataF[i];
      delete [] dataF;
      dataF = nullptr;
    }
    if (dataI) {
      for (i=0; i<DIMENSIONS; ++i) 
        if(dataI[i]) delete [] dataI[i];
      delete [] dataI;
      dataI = nullptr;
    }
  }

  void SimData::reserve(int num) {
    // Clear old pointers
    clean();

    // Set number and size
    number = 0;
    numberG = 0;
    size = num;

    // Reserve new arrays of arrays
    x = new RealType* [size];
    for (int n=0; n<size; ++n) x[n] = new RealType[DIMENSIONS];
    v = new RealType* [size];
    for (int n=0; n<size; ++n) v[n] = new RealType[DIMENSIONS];
    f = new RealType* [size];
    for (int n=0; n<size; ++n) f[n] = new RealType[DIMENSIONS];

    // Reserve new arrays
    sg   = new RealType[size];
    im   = new RealType[size];
    type = new int[size];

    // Set types to -1
    for (int n=0; n<size; ++n) type[n] = -1;
  }

  void SimData::reserveAll(int num) {
    // Clear old pointers
    clean();

    // Set number and size
    number = 0;
    numberG = 0;
    size = num;

    // Reserve new arrays of arrays
    x = new RealType* [size];
    for (int n=0; n<size; ++n) x[n] = new RealType[DIMENSIONS];
    v = new RealType* [size];
    for (int n=0; n<size; ++n) v[n] = new RealType[DIMENSIONS];
    f = new RealType* [size];
    for (int n=0; n<size; ++n) f[n] = new RealType[DIMENSIONS];

    // Reserve new arrays
    sg   = new RealType[size];
    im   = new RealType[size];
    type = new int[size];

    // Set types to -1
    for (int n=0; n<size; ++n) type[n] = -1;

    // Reserve for extra data
    dataF = new RealType*[size];
    for (int n=0; n<size; ++n) dataF[n] = nullptr;
    dataI = new int*[size];
    for (int n=0; n<size; ++n) dataI[n] = nullptr;
  }

  void SimData::clearX() {
    // Set all positions to zero
    for (int n=0; n<number; ++n) 
      for (int d=0; d<DIMENSIONS; ++d) x[n][d] = 0;
  }

  void SimData::clearV() {
    // Set all velocities to zero
    for (int n=0; n<number; ++n) 
      for (int d=0; d<DIMENSIONS; ++d) v[n][d] = 0;
  }

  void SimData::clearF() {
    // Set all forces to zero
    for (int n=0; n<number; ++n) 
      for (int d=0; d<DIMENSIONS; ++d) f[n][d] = 0;
  }
  
}