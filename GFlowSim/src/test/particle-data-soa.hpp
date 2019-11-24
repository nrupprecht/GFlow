#ifndef __PARTICLE_CONTAINER__GFLOW__
#define __PARTICLE_CONTAINER__GFLOW__

#include "../gflow.hpp"
#include "../other/timedobject.hpp"

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <unordered_map>

#include "../utility/memory.hpp"
#include "../utility/generic-dimension.hpp"
#include "../utility/vectormath.hpp"

#include "../utility/simd_generic.hpp"

#include "d-vec.hpp"
#include "container-base.hpp"

namespace GFlowSimulation {

  // template<int dims, DataLayout layout> class ParticleContainer : public ContainerBase {

  // };

  template<int dims> class ParticleContainer_SOA : public ContainerBase {
  public:
    //! \brief Constructor.
    ParticleContainer_SOA(GFlow*);

    //! \brief Destructor.
    ~ParticleContainer_SOA();

    //! \brief Initialize the particle container. After this point, no new entries should be added.
    virtual void initialize() override;

    //! \brief Resets timers.
    virtual void pre_integrate() override;

    //! \brief Remove all halo and ghost particles.
    virtual void post_integrate() override;

    //! \brief Is this the structure-of-arrays version (yes).
    bool is_soa();

    //! \brief Update simdata, migrate particles to other processors, handle assignment and initialization of ghost particles.
    void update();

    int add_vector_entry(const string& name);

    int add_scalar_entry(const string& name);

    int add_integer_entry(const string& name);

    // --- Accessing structs.
    struct vec_access;
    struct scalar_access;
    struct integer_access;

    vec_access X() { return vec_access(vector_data[0]); }
    vec_access V() { return vec_access(vector_data[1]); }
    vec_access F() { return vec_access(vector_data[2]); }
    vec_access v_entry(int i) { return vec_access(vector_data[i]); }

    vec<dims> X(int i) { return vec<dims>(vector_data[0][i]); }
    vec<dims> V(int i) { return vec<dims>(vector_data[1][i]); }
    vec<dims> F(int i) { return vec<dims>(vector_data[2][i]); }

    scalar_access R() { return scalar_access(scalar_data[0]); }
    scalar_access Im() { return scalar_access(scalar_data[1]); }
    scalar_access s_entry(int i) { return scalar_access(scalar_data[i]); }

    real& R(int i) { return scalar_data[0][i]; }
    real& Im(int i) { return scalar_data[1][i]; }

    integer_access Type() { return integer_access(integer_data[0]); }
    integer_access Id() { return integer_access(integer_data[1]); }
    integer_access i_entry(int i) { return integer_access(integer_data[i]); }

    int& Type(int i) { return integer_data[0][i]; }
    int& Id(int i) { return integer_data[1][i]; }

    //! \brief Clear all the [ar]-th vector entries of all the particles.
    void clear_vec(int ar);

    //! \brief Add an owned particle to the simulation (i.e. a particle that is owned by this processor) with zero velocity.
    //!
    //! Assumes that there are no ghost particles currently stored.
    int add_particle(vec<dims> x, real r, real mass, int type=0);

    //! \brief Add an owned particle to the simulation (i.e. a particle that is owned by this processor).
    //!
    //! Assumes that there are no ghost particles currently stored.
    int add_particle(vec<dims> x, vec<dims> v, real r, real mass, int type=0);

    //! \brief Reserve space for particles. 
    //!
    //! Assumes that there are no ghost particles. Only owned particles will be transfered (if there are any).
    void reserve(uint s);

    //! \brief Mark a particle for removal.
    void mark_for_removal(int id);

    //! \brief Remove particles that have been marked, and fill in space. Particles will contiguous after this function.
    void do_particle_removal();

  private:

    //! \brief Resize the particle data memory so that more (or fewer) particles fit in memory. Owned particles (up to total_size) are
    //! transfered to the new memory.
    void resize(const int total_size);

    //! \brief Vector data. Vectors are [dim] - dimensional vectors.
    vector<real**> vector_data;
    
    //! \brief Scalar data. Single entries.
    vector<real*> scalar_data;

    //! \brief Integer data.
    vector<int*> integer_data;
  };

  // Template function implementations.
  #include "particle-data-soa.tpp"
  #include "particle-data-soa-access.tpp"

};

#endif // __PARTICLE_CONTAINER__GFLOW__