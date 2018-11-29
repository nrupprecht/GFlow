#ifndef __ATOMS_HPP_GFLOW__
#define __ATOMS_HPP_GFLOW__

#include "../gflow.hpp"
#include "../utility/vectormath.hpp"

#include <set> // For storing sets of holes in the arrays
#include <map> // For mapping owned particles to halo particles

namespace GFlowSimulation {

  class Atoms : public Base {
  public:
    //! @brief Default constructor.
    Atoms(GFlow*);

    //! @brief Destructor.
    ~Atoms();

    //! @brief Initialize the atom container.
    virtual void initialize() override;

    //! @brief Reserve space for particles, extending the lengths of all arrays to the requested size.
    void reserve(int);

    //! @brief Add an owned particle to the simdata.
    //!
    //! Add a particle to the simdata. This is the public version of the function, so we can only add owned particles. 
    //! Depending on how many particles are in the array, and the array capacities, it may be necessary to resize the array 
    //! to add the particle.
    void addParticle(const RealType*, const RealType*, const RealType, const RealType, const int);

    //! @brief Exchange particles with neighboring domains.
    //void exchange_particles();

    // --- Accessors

    // --- Get vector data
    RealType** X();
    RealType*  X_arr();
    RealType*  X(int);
    RealType&  X(int, int);
    RealType** V();
    RealType*  V_arr();
    RealType*  V(int);
    RealType&  V(int, int);
    RealType** F();
    RealType*  F_arr();
    RealType*  F(int);
    RealType&  F(int, int);

    // --- Get scalar data
    RealType* Sg();
    RealType& Sg(int);
    RealType* Im();
    RealType& Im(int);

    // --- Get integer data
    int* Type();
    int& Type(int);
    int* Id();
    int& Id(int);

    //! @brief Return the number of owned particles.
    int size();
    int Number();

  private:
    // --- Helper functions

    //! @brief Allocate more space to hold owned particles.
    void resize_owned(int);

    //! @brief Vector data.
    //!
    //! Contains postion (0), velocity (1), and force (2).
    vector<RealType**> vdata;

    //! @brief Scalar data.
    //! 
    //! Contains sigma (0), inverse mass (1). Can also contain repulsion, dissipation, coefficient of friction, etc.
    vector<RealType*>  sdata;

    //! @brief Integer data.
    //!
    //! Contains type (0), global id (1). Can also contain body membership information, etc.
    vector<int*> idata;

    //! @brief The next global id a particle will be given.
    int next_global_id;

    //! @brief A map between global and local ids, <global, local>
    std::map<int, int> id_map;

    //! @brief Number of particles on this processor.
    int number; 

    //! @brief The index of the last spot that is reserved for owned particles.
    int last_owned;

    //! @brief The total number of entries in each nonnull array.
    int capacity;
  };

}


#endif // __ATOMS_HPP_GFLOW__