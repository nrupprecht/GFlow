#include "simdata.hpp"

namespace GFlowSimulation {

  SimData::SimData(GFlow *gflow) : Base(gflow), number(0), numberG(0), x(nullptr), v(nullptr), f(nullptr), 
                                  sg(nullptr), im(nullptr), type(nullptr), dataF(nullptr), dataI(nullptr)
  {
    for(int i=0; i<DIMENSIONS; ++i) wrap[i] = true;
  };

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
  
}