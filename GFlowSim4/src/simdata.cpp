#include "simdata.hpp"
// Other files
#include "memory.hpp"
#include "vectormath.hpp"
#include "memoryoptimizer.hpp"

namespace GFlowSimulation {

  SimData::SimData(GFlow *gflow) : Base(gflow), number(0), end_owned(0), number_halo(0), end_halo(0), 
    number_ghost(0), end_ghost(0), x(nullptr), v(nullptr), f(nullptr), sg(nullptr), im(nullptr), type(nullptr), 
    dataF(nullptr), dataI(nullptr), body(nullptr) {};

  SimData::~SimData() {
    clean();
  }

  void SimData::initialize() {
    // Call base initialization
    Base::initialize();
    // Reorder memory
    // MemoryOptimizer::GridParticles(*this, Base::gflow->getBounds());
  }

  void SimData::clean() {
    // Theses are arrays of arrays
    if (x) dealloc_array_2d(x);
    if (v) dealloc_array_2d(v);
    if (f) dealloc_array_2d(f);
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
      for (int i=0; i<DIMENSIONS; ++i) 
        if(dataF[i]) delete [] dataF[i];
      delete [] dataF;
      dataF = nullptr;
    }
    if (dataI) {
      for (int i=0; i<DIMENSIONS; ++i) 
        if(dataI[i]) delete [] dataI[i];
      delete [] dataI;
      dataI = nullptr;
    }
  }

  void SimData::reserve(int num) {
    // Clear old pointers
    clean();
    // At this point, dataF, dataI are both nullptr.

    // Don't allocate for a zero or negative number of particles
    if (num<=0) return;

    // Set numbers
    number = 0;
    number_halo = 0;
    number_ghost = 0;
    // ONLY reserving space for owned particles
    end_owned = num;
    end_halo = num;
    end_ghost = num;

    // Reserve new arrays of arrays
    x = alloc_array_2d<RealType>(num, DIMENSIONS);
    v = alloc_array_2d<RealType>(num, DIMENSIONS);
    f = alloc_array_2d<RealType>(num, DIMENSIONS);

    // Reserve new arrays
    sg   = new RealType[num];
    im   = new RealType[num];
    type = new int[num];

    // Set types to -1
    for (int n=0; n<num; ++n) type[n] = -1;
  }
  
  //! @param num The number of owned particles to reserve data space for.
  void SimData::reserveAll(int num) {
    // Clear old pointers
    clean();

    // Set numbers
    number = 0;
    number_halo = 0;
    number_ghost = 0;
    // ONLY reserving space for owned particles
    end_owned = num;
    end_halo = num;
    end_ghost = num;

    // Reserve new arrays of arrays
    x = alloc_array_2d<RealType>(num, DIMENSIONS);
    v = alloc_array_2d<RealType>(num, DIMENSIONS);
    f = alloc_array_2d<RealType>(num, DIMENSIONS);

    // Reserve new arrays
    sg   = new RealType[num];
    im   = new RealType[num];
    type = new int[num];

    // Set types to -1
    for (int n=0; n<num; ++n) type[n] = -1;

    // Reserve for extra data
    dataF = new RealType*[num];
    for (int n=0; n<num; ++n) dataF[n] = nullptr;
    dataI = new int*[num];
    for (int n=0; n<num; ++n) dataI[n] = nullptr;
  }

  //! @param num_owned The number of owned particles to make space for.
  //! @param num_halo  The number of halo particles to make space for.
  //! @param num_ghost The number of ghost particles to make space for.
  //! @param aux_F If true, we reserve space in dataF.
  //! @param aux_I If true, we reserve space in dataI.
  void SimData::reserve(int num_owned, int num_halo, int num_ghost, bool aux_F, bool aux_I) {

  }

  //! @param resize_owned The amount of extra space to reserve for owned particls.
  //! @param resize_halo  The amount of extra space to reserve for halo particles.
  //! @param resize_ghost The amount of extra space to reserve for ghost particles.
  void SimData::resize(int resize_owned, int resize_halo, int resize_ghost) {
    // Find the new total amount of space we need. End ghost is the size of the entire array
    int size = end_ghost + resize_owned + resize_halo + resize_ghost; 
    // Allocate new arrays
    RealType **nx  = alloc_array_2d<RealType>(size, DIMENSIONS);
    RealType **nv  = alloc_array_2d<RealType>(size, DIMENSIONS);
    RealType **nf  = alloc_array_2d<RealType>(size, DIMENSIONS);
    RealType *nsg = new RealType[size];
    RealType *nim  = new RealType[size];
    int *ntype = new int[size];
    RealType **ndataF = nullptr;
    int **ndataI = nullptr;
    if (dataF) ndataF = new RealType*[size];
    if (dataI) ndataI = new int*[size];
    // Copy over data
    copyHelper<DIMENSIONS>(resize_owned, resize_halo, resize_ghost, x[0], nx[0]);
    copyHelper<DIMENSIONS>(resize_owned, resize_halo, resize_ghost, v[0], nv[0]);
    copyHelper<DIMENSIONS>(resize_owned, resize_halo, resize_ghost, f[0], nf[0]);
    copyHelper<1>(resize_owned, resize_halo, resize_ghost, sg, nsg);
    copyHelper<1>(resize_owned, resize_halo, resize_ghost, im, nim);
    copyHelper<1>(resize_owned, resize_halo, resize_ghost, type, ntype);
    // Potentially copy auxilary data
    if (dataF) copyHelper<1>(resize_owned, resize_halo, resize_ghost, dataF, ndataF);
    if (dataF) copyHelper<1>(resize_owned, resize_halo, resize_ghost, dataI, ndataI);
  }

  void SimData::clearX() {
    // Set all positions to zero
    RealType *_x = x[0];
    for (int i=0; i<number*DIMENSIONS; ++i) _x[i] = 0;
  }

  void SimData::clearV() {
    // Set all velocities to zero
    RealType *_v = v[0];
    for (int i=0; i<number*DIMENSIONS; ++i) _v[i] = 0;
  }

  void SimData::clearF() {
    // Set all forces to zero
    RealType *_f = f[0];
    for (int i=0; i<number*DIMENSIONS; ++i) _f[i] = 0;
  }

  //! @param x The position of the particle.
  //! @param v The velocity of the particle.
  //! @param sg The cutoff radius of the particle.
  //! @param im The inverse mass of the particle.
  //! @param type The type of the particle.
  void SimData::addParticle(const RealType *X, const RealType *V, const RealType Sg, 
    const RealType Im, const int Type) 
  {
    if (end_owned<=number) 
      resize(ceil(0.25*number), 0, 0);
    // Copy data
    copyVec(X, x[number]);
    copyVec(V, v[number]);
    sg[number]   = Sg;
    im[number]   = Im;
    type[number] = Type;
    // There is now one more particle
    ++number;
    // Set flag
    needs_remake = true;
  }

  //! @param X The position of the particle.
  //! @param V The velocity of the particle.
  //! @param Sg The cutoff radius of the particle.
  //! @param Im The inverse mass of the particle.
  //! @param Type The type of the particle.
  //! @param own_type What type of ownership the processor has over this particle. @see ParticleType
  void SimData::addParticle(const RealType *X, const RealType *V, const RealType Sg, 
    const RealType Im, const int Type, const ParticleOwnership own_type) 
  {
    switch (own_type) {
      default:
      case ParticleOwnership::Owned: {
        // Check if we need to resize the array
        if (end_owned<=number) 
          resize(ceil(0.25*number), 0, 0);
        // Copy data
        copyVec(X, x[number]);
        copyVec(V, v[number]);
        sg[number]   = Sg;
        im[number]   = Im;
        type[number] = Type;
        // There is now one more particle
        ++number;
        // Done.
        break;
      }
      case ParticleOwnership::Halo: {
        // Check if we need to resize the array
        if (end_halo - end_owned <= number_halo) resize(0, ceil(0.25*number_halo), 0);
        // Copy data
        copyVec(X, x[number_halo+end_owned]);
        copyVec(V, v[number_halo+end_owned]);
        sg[number_halo+end_owned]   = Sg;
        im[number_halo+end_owned]   = Im;
        type[number_halo+end_owned] = Type;
        // There is now one more halo particle
        ++number_halo;
        // Done.
        break;
      }
      case ParticleOwnership::Ghost: {
        // Check if we need to resize the array
        if (end_ghost - end_halo <= number_ghost) resize(0, 0, ceil(0.25*number_ghost));
        // Copy data
        copyVec(X, x[number_ghost+end_halo]);
        copyVec(V, v[number_ghost+end_halo]);
        sg[number_ghost+end_halo]   = Sg;
        im[number_ghost+end_halo]   = Im;
        type[number_ghost+end_halo] = Type;
        // There is now one more ghost particle
        ++number_ghost;
        // Done.
        break;
      }
    }
    // Set flag
    needs_remake = true;
  }

  void SimData::addHaloParticle(const RealType *X, const RealType *V, const RealType Sg, 
    const RealType Im, const int Type)
  {
    // Check if we need to resize the array
    if (end_halo - end_owned <= number_halo) resize(0, ceil(0.25*number_halo), 0);
    // Copy data
    copyVec(X, x[number_halo+end_owned]);
    copyVec(V, v[number_halo+end_owned]);
    sg[number_halo+end_owned]   = Sg;
    im[number_halo+end_owned]   = Im;
    type[number_halo+end_owned] = Type;
    // There is now one more halo particle
    ++number_halo;
    // Set flag
    needs_remake = true;
  }

  void SimData::addGhostParticle(const RealType *X, const RealType *V, const RealType Sg, 
    const RealType Im, const int Type)
  {
    // Check if we need to resize the array
    if (end_ghost - end_halo <= number_ghost) resize(0, 0, ceil(0.25*number_ghost));
    // Copy data
    copyVec(X, x[number_ghost+end_halo]);
    copyVec(V, v[number_ghost+end_halo]);
    sg[number_ghost+end_halo]   = Sg;
    im[number_ghost+end_halo]   = Im;
    type[number_ghost+end_halo] = Type;
    // There is now one more ghost particle
    ++number_ghost;
    // Set flag
    needs_remake = true;
  }

  void SimData::clearHaloParticles() {
    for (int i=end_owned; i<end_owned+number_halo; ++i) type[i] = -1;
    number_halo = 0;
  }

  void SimData::clearGhostParticles() {
    for (int i=end_halo; i<end_halo+number_ghost; ++i) type[i] = -1;
    number_ghost = 0;
  }
  
}
