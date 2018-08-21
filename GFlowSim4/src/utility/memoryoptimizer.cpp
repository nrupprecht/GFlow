#include "memoryoptimizer.hpp"
// Other files
#include "../base/simdata.hpp"
#include "array.hpp"

namespace GFlowSimulation {

  void MemoryOptimizer::GridParticles(SimData& simData, const Bounds& bounds) {
    // Check that there's something to do 
    int number = simData.number;
    if (number<=0) return;

    // Compute sizes
    int sizes[DIMENSIONS];
    RealType width[DIMENSIONS];
    int total = 1;
    for (int d=0; d<DIMENSIONS; ++d) {
      sizes[d] = 4;
      width[d] = (bounds.max[d] - bounds.min[d])/sizes[d];
      total *= sizes[d];
    }

    // An array we can use to sort particles
    Array<vector<int> > grid(sizes);

    auto getLinearIndex = [&] (RealType *x)->int {
      int index = 0;
      int product = 1;
      for (int d=0; d<DIMENSIONS; ++d) {
        int X = static_cast<int>((x[d] - bounds.min[d])/width[d]);
      	X = X<0 ? 0 : X;
      	X = X>=sizes[d] ? sizes[d]-1 : X;
        index += (X*product); 
        product *= sizes[d];
      }
      return index;
    };

    // Sort particles into boxes

    for (int i=0; i<number; ++i) {
      int index = getLinearIndex(simData.x[i]);
      grid[index].push_back(i);
    }

    // --- Reorder them in their lists box by box.
    // Reserve new arrays and arrays-of-arrays
    RealType **x = alloc_array_2d<RealType>(number, DIMENSIONS);
    RealType **v = alloc_array_2d<RealType>(number, DIMENSIONS);
    RealType **f = alloc_array_2d<RealType>(number, DIMENSIONS);
    RealType *sg   = new RealType[number];
    RealType *im   = new RealType[number];
    int      *type = new int[number];
    //! @todo Copy auxilary data too.

    // Sort
    int count = 0;
    for (int i=0; i<total; ++i) {
      for (auto id : grid[i]) {
        copyParticle(simData, id, x[count], v[count], f[count], sg[count], im[count], type[count]);
        ++count;
      }
    }

    // --- Destroy old data, set new data
    // Position
    dealloc_array_2d(simData.x);
    simData.x = x;
    // Velcity
    dealloc_array_2d(simData.v);
    simData.v = v;
    // Force
    dealloc_array_2d(simData.f);
    simData.f = f;
    // Radius
    delete [] simData.sg;
    simData.sg = sg;
    // Inverse mass
    delete [] simData.im;
    simData.im = im;
    // Particle type
    delete [] simData.type;
    simData.type = type;
  }

}
