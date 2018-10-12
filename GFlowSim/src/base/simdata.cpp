#include "simdata.hpp"
// Other files
#include "../utility/memory.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  SimData::SimData(GFlow *gflow) : Base(gflow) {
    data_array = vector<RealType*>(3, nullptr);
  };

  SimData::~SimData() {
    clean();
  }

  void SimData::initialize() {
    // Call base initialization
    Base::initialize();
  }

  void SimData::clean() {
    // Theses are arrays of arrays
    if (x) dealloc_array_2d(x);
    if (v) dealloc_array_2d(v);
    if (f) dealloc_array_2d(f);
    clearAngularDynamics();
    // Clear data arrays
    for (auto &dt : data_array)
      if (dt) delete [] dt;
    if (type) {
      delete [] type;
      type = nullptr;
    }
    // --- Delete additional data - these are vectors of arrays
    for (auto &ar : dataF)
      if (ar) {
        delete [] ar;
        ar = nullptr;
      }
    for (auto &ar : dataI)
      if (ar) {
        delete [] ar;
        ar = nullptr;
      }
  }

  void SimData::reserve(int num) {
    // Clear old pointers
    clean();
    // At this point, dataF, dataI are both empty.

    // Don't allocate for a zero or negative number of particles
    if (num<=0) return;

    // Set numbers
    number       = 0;
    number_halo  = 0;
    number_ghost = 0;
    // ONLY reserving space for owned particles
    end_owned = num;
    end_halo = num;
    end_ghost = num;

    // Reserve new arrays of arrays
    x = alloc_array_2d<RealType>(num, DIMENSIONS);
    v = alloc_array_2d<RealType>(num, DIMENSIONS);
    f = alloc_array_2d<RealType>(num, DIMENSIONS);
    if (angularDynamics) allocateAngularDynamics(num);

    // Reserve new arrays
    /*
    sg   = new RealType[num];
    im   = new RealType[num];
    */
    data_array[0] = new RealType[num]; // Sigmas
    data_array[1] = new RealType[num]; // Inverse masses
    type = new int[num];
    global_id = new int[num];

    // Set types to -1
    for (int n=0; n<num; ++n) type[n] = -1;
  }
  
  //! @param num The number of owned particles to reserve data space for.
  void SimData::reserveAll(int num) {
    // Clear old pointers
    clean();

    // Set numbers
    number       = 0;
    number_halo  = 0;
    number_ghost = 0;
    // ONLY reserving space for owned particles
    end_owned = num;
    end_halo = num;
    end_ghost = num;

    // Reserve new arrays of arrays
    x = alloc_array_2d<RealType>(num, DIMENSIONS);
    v = alloc_array_2d<RealType>(num, DIMENSIONS);
    f = alloc_array_2d<RealType>(num, DIMENSIONS);

    if (angularDynamics) allocateAngularDynamics(num);

    // Reserve new arrays
    /*
    sg   = new RealType[num];
    im   = new RealType[num];
    */
    data_array[0] = new RealType[num]; // Sigmas
    data_array[1] = new RealType[num]; // Inverse masses
    type = new int[num];
    global_id = new int[num];

    // Set types to -1
    for (int n=0; n<num; ++n) type[n] = -1;

    // Reserve for extra data
    for (int i=0; i<dataF.size(); ++i)
      dataF[i] = new RealType[num];
    for (int i=0; i<dataI.size(); ++i)
      dataI[i] = new int[num];
  }

  //! @param num_owned The number of owned particles to make space for.
  //! @param num_halo  The number of halo particles to make space for.
  //! @param num_ghost The number of ghost particles to make space for.
  //! @param aux_F If true, we reserve space in dataF.
  //! @param aux_I If true, we reserve space in dataI.
  void SimData::reserve(int num_owned, int num_halo, int num_ghost, bool aux_F, bool aux_I) {
    throw Unimplemented();
  }

  //! @param resize_owned The amount of extra space to reserve for owned particls.
  //! @param resize_halo  The amount of extra space to reserve for halo particles.
  //! @param resize_ghost The amount of extra space to reserve for ghost particles.
  void SimData::resize(int resize_owned, int resize_halo, int resize_ghost) {
    // Find the new total amount of space we need. End ghost is (currently) the size of the entire array.
    int size = end_ghost + resize_owned + resize_halo + resize_ghost; 
    // Allocate new arrays
    RealType **nx  = alloc_array_2d<RealType>(size, DIMENSIONS);
    RealType **nv  = alloc_array_2d<RealType>(size, DIMENSIONS);
    RealType **nf  = alloc_array_2d<RealType>(size, DIMENSIONS);
    // Angular arrays
    RealType *nth = nullptr, *nom = nullptr, *ntq = nullptr, *niI = nullptr;
    if (angularDynamics && DIMENSIONS>1) {
      nth = new RealType[size];
      nom = new RealType[size];
      ntq = new RealType[size];
      niI = new RealType[size];
    }    
    RealType *nsg  = new RealType[size];
    RealType *nim  = new RealType[size];
    int *ntype     = new int[size];
    int *nbody     = new int[size];
    int *ngid      = new int[size];

    vector<RealType*> ndataF;
    vector<int*> ndataI;

    if (!dataF.empty()) {
      ndataF = vector<RealType*>(dataF.size(), nullptr);
      for (auto &ar : ndataF) ar = new RealType[end_ghost];
    }
    if (!dataI.empty()) {
      ndataI = vector<int*>(dataI.size(), nullptr);
      for (auto &ar : ndataI) ar = new int[end_ghost];
    }
    // --- Copy over data
    if (x) copyHelper<DIMENSIONS>(resize_owned, resize_halo, resize_ghost, x[0], nx[0]);
    if (v) copyHelper<DIMENSIONS>(resize_owned, resize_halo, resize_ghost, v[0], nv[0]);
    if (f) copyHelper<DIMENSIONS>(resize_owned, resize_halo, resize_ghost, f[0], nf[0]);
    if (data_array[0]) copyHelper<1>(resize_owned, resize_halo, resize_ghost, data_array[0], nsg);
    if (data_array[1]) copyHelper<1>(resize_owned, resize_halo, resize_ghost, data_array[1], nim);
    if (type) copyHelper<1>(resize_owned, resize_halo, resize_ghost, type, ntype);
    if (body) copyHelper<1>(resize_owned, resize_halo, resize_ghost, body, nbody);
    if (global_id) copyHelper<1>(resize_owned, resize_halo, resize_ghost, global_id, ngid);
    // --- Copy over angular data
    if (angularDynamics && DIMENSIONS>1) {
      copyHelper<1>(resize_owned, resize_halo, resize_ghost, th, nth);
      copyHelper<1>(resize_owned, resize_halo, resize_ghost, om, nom);
      copyHelper<1>(resize_owned, resize_halo, resize_ghost, tq, ntq);
      copyHelper<1>(resize_owned, resize_halo, resize_ghost, data_array[2], niI);
    }
    // --- Copy auxilary data
    for (int i=0; i<dataF.size(); ++i)
      copyHelper<1>(resize_owned, resize_halo, resize_ghost, dataF[i], ndataF[i]);
    for (int i=0; i<dataI.size(); ++i)
      copyHelper<1>(resize_owned, resize_halo, resize_ghost, dataI[i], ndataI[i]);
    // --- Delete old arrays, set new ones
    clean();
    // Old arrays are cleaned, safe to set the new ones
    x  = nx;
    v  = nv;
    f  = nf;
    th = nth;
    om = nom;
    tq = ntq;
    data_array[0] = nsg;
    data_array[1] = nim;
    data_array[2] = niI;
    type  = ntype;
    body  = nbody;
    dataF = ndataF;
    dataI = ndataI;
    global_id = ngid;
    // Update end_XYZ variables
    end_owned += resize_owned;
    end_halo  += (resize_owned + resize_halo);
    end_ghost += (resize_owned + resize_halo + resize_ghost);
  }

  void SimData::clearX() {
    // Check if there are particles
    if (number==0) return;
    // Set all positions to zero
    RealType *_x = x[0];
    for (int i=0; i<number*DIMENSIONS; ++i) _x[i] = 0.;
  }

  void SimData::clearV() {
    // Check if there are particles
    if (number==0) return;
    // Set all velocities to zero
    RealType *_v = v[0];
    for (int i=0; i<number*DIMENSIONS; ++i) _v[i] = 0.;
  }

  void SimData::clearF() {
    // Check if there are particles
    if (number==0) return;
    // Set all forces to zero
    RealType *_f = f[0];
    for (int i=0; i<number*DIMENSIONS; ++i) _f[i] = 0.;
  }

  void SimData::clearTh() {
    // Check if there are particles and angular variables are allocated
    if (number==0 || !angularDynamics || DIMENSIONS==1) return;
    // Set all forces to zero
    for (int i=0; i<number; ++i) th[i] = 0.;
  }

  void SimData::clearOm() {
    // Check if there are particles and angular variables are allocated
    if (number==0 || !angularDynamics || DIMENSIONS==1) return;
    // Set all omegas to zero
    for (int i=0; i<number; ++i) om[i] = 0.;
  }

  void SimData::clearTq() {
    // Check if there are particles and angular variables are allocated
    if (number==0 || !angularDynamics || DIMENSIONS==1) return;
    // Set all torques to zero
    for (int i=0; i<number; ++i) tq[i] = 0.;
  }

  void SimData::cprojF(RealType* dir) {
    // Make sure dir is normalized
    RealType norm[DIMENSIONS];
    normalVec(dir, norm);

    // Check if the direction is a standard coordinate direction
    for (int d=0; d<DIMENSIONS; ++d) {
      if (norm[d]==1.) {

        RealType *_f = f[0];
        for (int i=d; i<number*DIMENSIONS; i+=DIMENSIONS) _f[i] = 0.;
        return;
      }
    }

    // If the direction is arbitrary
    RealType mag = 0;
    RealType fproj[DIMENSIONS];
    for (int i=0; i<number; ++i) {
      mag = dotVec(f[i], norm);
      scalarMultVec(mag, f[i], fproj);
      minusEqVec(f[i], fproj);
    }
  }

  void SimData::setAllDataF(string entry_name, RealType value) {
    int id = getDataFEntry(entry_name);
    if (id>-1) {
      for (int i=0; i<number; ++i)
        dataF[id][i] = value;
    }
  }

  void SimData::setAllDataI(string entry_name, int value) {
    int id = getDataIEntry(entry_name);
    if (id>-1) {
      for (int i=0; i<number; ++i)
        dataI[id][i] = value;
    }
  }

  const BCFlag* SimData::getBCs() const {
    return Base::gflow->getBCs();
  }

  int SimData::getLocalID(int global) const {
    auto it = id_map.find(global);
    // Return the global iterator. We use -1 to mean "no such particle."
    return it==id_map.end() ? -1 : it->second;
  }

  int SimData::getNextGlobalID() const {
    return next_global_id;
  }

  Bounds SimData::getBounds() const {
    return Base::gflow->getBounds();
  }

  bool SimData::usingAngularDynamics() const {
    return angularDynamics;
  }

  //! @param x The position of the particle.
  //! @param v The velocity of the particle.
  //! @param sg The cutoff radius of the particle.
  //! @param im The inverse mass of the particle.
  //! @param type The type of the particle.
  void SimData::addParticle(const RealType *X, const RealType *V, const RealType Sg, const RealType Im, const int Type) 
  {
    if (end_owned<=number) {
      // Calculate the new amount of size we want
      int resize_owned = max(static_cast<int>(ceil(0.25*number)), 32);
      resize(resize_owned, 0, 0);
    }
    // Copy data
    copyVec(X, x[number]);
    copyVec(V, v[number]);
    zeroVec(f[number]);
    // Insert into the global id map
    id_map.insert(IPair(next_global_id, number));
    // Assign a global id - this assumes that this is a *new* particle
    global_id[number] = ++next_global_id;
    data_array[0][number]   = Sg;
    data_array[1][number]   = Im;
    type[number] = Type;
    // If using angular dynamics, but iI is unspecified, treat as uniform disk.
    if (angularDynamics)
      data_array[2][number] = 0.5*(1./Im)*sqr(Sg);
    // There is now one more particle
    ++number;
    // Set flag
    needs_remake = true;
  }

  void SimData::addParticle(const RealType *X, const RealType *V, const RealType Sg, const RealType Im, const int Type, const RealType II) {
    if (!angularDynamics) addParticle(X, V, Sg, Im, Type);
    else {
      if (end_owned<=number) {
        // Calculate the new amount of size we want
        int resize_owned = max(static_cast<int>(ceil(0.25*number)), 32);
        resize(resize_owned, 0, 0);
      }
      // Copy data
      copyVec(X, x[number]);
      copyVec(V, v[number]);
      zeroVec(f[number]);
      // Insert into the global id map
      id_map.insert(IPair(next_global_id, number));
      // Assign a global id - this assumes that this is a *new* particle
      global_id[number] = ++next_global_id;
      data_array[0][number]   = Sg;
      data_array[1][number]   = Im;
      data_array[2][number]   = II;
      type[number] = Type;
      // There is now one more particle
      ++number;
      // Set flag
      needs_remake = true;
    }
  }

  //! @param X The position of the particle.
  //! @param V The velocity of the particle.
  //! @param Sg The cutoff radius of the particle.
  //! @param Im The inverse mass of the particle.
  //! @param Type The type of the particle.
  //! @param own_type What type of ownership the processor has over this particle. @see ParticleType
  void SimData::addParticle(const RealType *X, const RealType *V, const RealType Sg, const RealType Im, const int Type, const ParticleOwnership own_type) 
  {
    // @todo Global ID transfer
    switch (own_type) {
      default:
      case ParticleOwnership::Owned: {
        // Use the normal add particle function
        addParticle(X, V, Sg, Im, Type);
        // Done.
        break;
      }
      case ParticleOwnership::Halo: {
      	// Temporary
      	throw false;
        // Check if we need to resize the array
        if (end_halo - end_owned <= number_halo) resize(0, ceil(0.25*number_halo), 0);
        // Copy data
        copyVec(X, x[number_halo+end_owned]);
        copyVec(V, v[number_halo+end_owned]);
        data_array[0][number_halo+end_owned]   = Sg;
        data_array[1][number_halo+end_owned]   = Im;
        type[number_halo+end_owned] = Type;
        // There is now one more halo particle
        ++number_halo;
        // Done.
        break;
      }
      case ParticleOwnership::Ghost: {
      	// Temporary
      	throw false;
        // Check if we need to resize the array
        if (end_ghost - end_halo <= number_ghost) resize(0, 0, ceil(0.25*number_ghost));
        // Copy data
        copyVec(X, x[number_ghost+end_halo]);
        copyVec(V, v[number_ghost+end_halo]);
        data_array[0][number_ghost+end_halo]   = Sg;
        data_array[1][number_ghost+end_halo]   = Im;
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

  void SimData::markForRemoval(const int id) {
    remove_list.insert(id);
  }

  void SimData::doParticleRemoval() {
    // If there is nothing to remove, we're done
    if (remove_list.empty() || number==0) return;
    // Variables
    int count = 0, count_back = number, need_removal = 0;
    // Set all types to -1, remove global ids
    for (auto id : remove_list) {
      if (type[id]!=-1) ++need_removal;
      type[id] = -1;
    }
    // Fill in all holes
    int removed = 0;
    for(auto id : remove_list) {
      // We either need to start at (number-1), or moved the particle that was at count_back. Either way, decrement.
      --count_back;
      // We have removed a particle
      ++removed;

      // Find the next valid particle (counting back from the end) to fill for the particle we want to remove
      // C++ 20 has a std::set contains() function.
      while ( contains(remove_list, count_back) && count_back>id) --count_back;

      // Move the particle to fill the particle we want to remove - moving the good particle to the hole erases the hole's 
      // global id.
      if (count_back>id) moveParticle(count_back, id);
      else break;
    }
    // Decrease number
    number -= need_removal;

    // Clear list
    remove_list.clear();
    // We need to update
    if (removed>0) needs_remake = true;
  }

  //! @param id The id (place in the data lists) of the particle that should be removed.
  void SimData::removeParticle(int id) {
    // We assume that this happens infrequently, so fill in the hole now. Move the last particle
    // in the array to this place. We must remake the verlet list though.
    if (id>=number || id<0) return; // Not a valid spot
    if (number>1) {
      moveParticle(number-1, id); // Moving the particle erases the removed particle's global id
      // We need to remake
      needs_remake = true;
    }
    else { // This was the only particle
      type[0] = -1;
      // We don't actually need to remake
    }
  }

  void SimData::addHaloParticle(const RealType *X, const RealType *V, const RealType Sg, 
    const RealType Im, const int Type)
  {
    // @todo Global ID transfer

    // Check if we need to resize the array
    if (end_halo - end_owned <= number_halo) resize(0, ceil(0.25*number_halo), 0);
    // Copy data
    copyVec(X, x[number_halo+end_owned]);
    copyVec(V, v[number_halo+end_owned]);
    data_array[0][number_halo+end_owned]   = Sg;
    data_array[1][number_halo+end_owned]   = Im;
    type[number_halo+end_owned] = Type;
    // There is now one more halo particle
    ++number_halo;
    // Set flag
    needs_remake = true;
  }

  void SimData::addGhostParticle(const RealType *X, const RealType *V, const RealType Sg, 
    const RealType Im, const int Type)
  {
    // @todo Global ID transfer

    // Check if we need to resize the array
    if (end_ghost - end_halo <= number_ghost) resize(0, 0, ceil(0.25*number_ghost));
    // Copy data
    copyVec(X, x[number_ghost+end_halo]);
    copyVec(V, v[number_ghost+end_halo]);
    data_array[0][number_ghost+end_halo]   = Sg;
    data_array[1][number_ghost+end_halo]   = Im;
    type[number_ghost+end_halo] = Type;
    // There is now one more ghost particle
    ++number_ghost;
    // Set flag
    needs_remake = true;
  }

  template<int width, typename T> inline void SimData::copyHelper(int resize_owned, int resize_halo, int resize_ghost, T *old_array, T *new_array) {
    // Move owned particle data
    for (int i=0; i<number*width; ++i) 
      new_array[i] = old_array[i];
    // Move halo data
    for (int i=end_owned*width; i<(end_owned+number_halo)*width; ++i) 
      new_array[i+resize_owned*width] = old_array[i];
    // Move ghost data
    for (int i=end_halo*width; i<(end_halo+number_ghost)*width; ++i)
      new_array[i+(resize_owned+resize_halo)*width] = old_array[i];
  }

  void SimData::moveParticle(int id_source, int id_dest) {
    // Get global ids
    int gs = global_id[id_source];
    int gd = global_id[id_dest];

    // --- Copy main data
    copyVec(x[id_source], x[id_dest]);
    copyVec(v[id_source], v[id_dest]);
    copyVec(f[id_source], f[id_dest]);
    global_id[id_dest] = gs;
    data_array[0][id_dest] = data_array[0][id_source];
    data_array[1][id_dest] = data_array[1][id_source];
    type[id_dest] = type[id_source];
    // --- Copy angular data
    if (angularDynamics) {
      th[id_dest] = th[id_source];
      om[id_dest] = om[id_source];
      tq[id_dest] = tq[id_source];
      data_array[2][id_dest] = data_array[2][id_source];
    }
    // --- Move auxilary data
    for (int i=0; i<dataF.size(); ++i)
      dataF[i][id_dest] = dataF[i][id_source];
    for (int i=0; i<dataI.size(); ++i)
      dataI[i][id_dest] = dataI[i][id_source];
    
    // Update id map - erases global id of the overwritten particle
    auto it = id_map.find(gs);
    if (id_map.end()!=it) it->second = id_dest;
    id_map.erase(gd);

    // Set id of source to -1
    type[id_source] = -1;
  }

  void SimData::clearHaloParticles() {
    for (int i=end_owned; i<end_owned+number_halo; ++i) type[i] = -1;
    number_halo = 0;
  }

  void SimData::clearGhostParticles() {
    for (int i=end_halo; i<end_halo+number_ghost; ++i) type[i] = -1;
    number_ghost = 0;
  }

  //! @brief Returns whether the simdata needs to be remade or not
  bool SimData::getNeedsRemake() { 
    return needs_remake; 
  }

  //! @brief Set the needs_remake flag
  void SimData::setNeedsRemake(bool r) { 
    needs_remake = r; 
  }

  void SimData::setAngularDynamics(bool ad) {
    if (angularDynamics && !ad) clearAngularDynamics();
    else if (!angularDynamics && ad) allocateAngularDynamics(number);
    angularDynamics = ad;
  }
  
  void SimData::addDataFEntry(const string entry_name, float default_value) {
    // Allocate the new array
    RealType *new_array = nullptr;
    if (end_ghost>0) new_array = new RealType[end_ghost];
    dataF.push_back(new_array);
    // Find the place of the array we just added
    int place = dataF.size()-1;
    // Set default value
    for (int i=0; i<end_ghost; ++i) dataF[place][i] = default_value;
    // Store name
    dataFIndex.insert(pair<string, int>(entry_name, place));
  }

  void SimData::addDataIEntry(const string entry_name, int default_value) {
    // Allocate the new array
    int *new_array = nullptr;
    if (end_ghost>0) new_array = new int[end_ghost];
    dataI.push_back(new_array);
    // Find the place of the array we just added
    int place = dataF.size()-1;
    // Set default value
    for (int i=0; i<end_ghost; ++i) dataI[place][i] = default_value;
    // Store name
    dataFIndex.insert(pair<string, int>(entry_name, place));
  }

  int SimData::getDataFEntry(const string entry_name) {
    auto p = dataFIndex.find(entry_name);
    if (p!=dataFIndex.end()) return p->second;
    return -1;
  }

  int SimData::getDataIEntry(const string entry_name) {
    auto p = dataIIndex.find(entry_name);
    if (p!=dataIIndex.end()) return p->second;
    return -1;
  }

  void SimData::clearAngularDynamics() {
    if (th) delete [] th;
    th = nullptr;
    if (om) delete [] om;
    om = nullptr;
    if (tq) delete [] tq;
    tq = nullptr;
    if (data_array[2]) delete [] data_array[2];
    data_array[2] = nullptr;
  }

  void SimData::allocateAngularDynamics(int num) {
    if (num<=0) return;
    // We only need one angle in 2d
    th = new RealType[num];
    om = new RealType[num];
    tq = new RealType[num];
    data_array[2] = new RealType[num];
  }

}
