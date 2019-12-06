#ifndef __PARTICLE_DATA_AOS_HPP__GFLOW__
#define __PARTICLE_DATA_AOS_HPP__GFLOW__

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

  template<int dims> class ParticleContainer<dims, DataLayout::AOS> : public ContainerBase {
  public:

    typedef ParticleContainer<dims, DataLayout::AOS> SelfType;

    //! \brief Default constructor.
    ParticleContainer(GFlow*);

    //! \brief Destructor.
    ~ParticleContainer();

    //! \brief Initialize the particle container. After this point, no new entries should be added.
    virtual void initialize() override;

    //! \brief Resets timers.
    virtual void pre_integrate() override;

    //! \brief Remove all halo and ghost particles.
    virtual void post_integrate() override;

    //! \brief Is this the structure-of-arrays version (no).
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

    vec_access X() { return vec_access(data_ptr, data_width, 0); }
    vec_access V() { return vec_access(data_ptr, data_width, dims); }
    vec_access F() { return vec_access(data_ptr, data_width, 2*dims); }
    vec_access v_entry(int i) { return vec_access(data_ptr, data_width, i*dims); }

    vec<dims> X(int i) { return vec<dims>(&data_ptr[i*data_width]); }
    vec<dims> V(int i) { return vec<dims>(&data_ptr[i*data_width + dims]); }
    vec<dims> F(int i) { return vec<dims>(&data_ptr[i*data_width + 2*dims]); }

    scalar_access R() { return scalar_access(data_ptr, data_width, n_vectors*dims); }
    scalar_access Im() { return scalar_access(data_ptr, data_width, n_vectors*dims + 1); }
    scalar_access s_entry(int i) { return scalar_access(data_ptr, data_width, n_vectors*dims + i); }

    real& R(int i) { return data_ptr[data_width*i + n_vectors*dims]; }
    real& Im(int i) { return data_ptr[data_width*i + n_vectors*dims + 1]; }

    integer_access Type() { return integer_access(data_ptr, data_width, n_vectors*dims + n_scalars); }
    integer_access Id() { return integer_access(data_ptr, data_width, n_vectors*dims + n_scalars + 1); }
    integer_access i_entry(int i) { return integer_access(data_ptr, data_width, n_vectors*dims + n_scalars + i); }

    void setType(int i, int ty) { data_ptr[i*data_width + n_vectors*dims + n_scalars] = *reinterpret_cast<real*>(&ty); }
    const int Type(int i) { return *reinterpret_cast<int*>(&data_ptr[i*data_width + n_vectors*dims + n_scalars]); }
    void setId(int i, int id) { data_ptr[i*data_width + n_vectors*dims + n_scalars + 1] = *reinterpret_cast<real*>(&id); }
    const int Id(int i) { return *reinterpret_cast<int*>(data_ptr[i*data_width + n_vectors*dims + n_scalars + 1]); }

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

  //private:

    //! \brief Resize the particle data memory so that more (or fewer) particles fit in memory. Owned particles (up to total_size) are
    //! transfered to the new memory.
    void resize(const int total_size);

    //! \brief The data width (in reals) of a single particle.
    int data_width = 3*dims + 4;

    //! \brief All data is stored in one array, with a particle made up of contiguous data entries.
    real *data_ptr = nullptr;
  };

  // Template function implementations.
  #include "particle-data-aos.tpp"
  #include "particle-data-aos-access.tpp"

}
#endif // __PARTICLE_DATA_AOS_HPP__GFLOW__