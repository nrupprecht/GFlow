#ifndef __PARTICLE_CONTAINER_BASE_HPP__GFLOW__
#define __PARTICLE_CONTAINER_BASE_HPP__GFLOW__

namespace GFlowSimulation {

  template<int dims> class ParticleContainer {

    virtual void initialize() = 0;

    int add_vector_entry(const string& name) = 0;
    int add_scalar_entry(const string& name) = 0;
    int add_integer_entry(const string& name) = 0;

    vec_access<dims> X() = 0;
    vec_access<dims> V() = 0;
    vec_access<dims> F() = 0;
    vec_access<dims> v_entry(int i) = 0;

    vec<dims> X(int i) = 0;
    vec<dims> V(int i) = 0;
    vec<dims> F(int i) = 0;

    scalar_access<dims> R() = 0;
    scalar_access<dims> Im() = 0;
    scalar_access<dims> s_entry(int i) = 0;

    real& R(int i) = 0;
    real& Im(int i) = 0;

    integer_access<dims> Type() = 0;
    integer_access<dims> Id() = 0;
    integer_access<dims> i_entry(int i) = 0;

    int& Type(int i) = 0;
    int& Id(int i) = 0;

    //! \brief Clear all the [ar]-th vector entries of all the particles.
    void clear_vec(int ar) = 0;

    //! \brief Add an owned particle to the simulation (i.e. a particle that is owned by this processor) with zero velocity.
    //!
    //! Assumes that there are no ghost particles currently stored.
    int add_particle(vec<dims> x, real r, real mass, int type=0) = 0;

    //! \brief Add an owned particle to the simulation (i.e. a particle that is owned by this processor).
    //!
    //! Assumes that there are no ghost particles currently stored.
    int add_particle(vec<dims> x, vec<dims> v, real r, real mass, int type=0) = 0;

    //! \brief Reserve space for particles. 
    //!
    //! Assumes that there are no ghost particles. Only owned particles will be transfered (if there are any).
    void reserve(uint s) = 0;

    //! \brief The size particle entries that may contain valid owned particles.
    int size() = 0;

    //! \brief The number of (valid) owned particles stored in this data structure.
    int number() = 0;

    void mark_for_removal(int id) = 0;

    //! \brief Remove particles that have been marked, and fill in space. Particles will contiguous after this function.
    void do_particle_removal() = 0;

  };
}
#endif // __PARTICLE_CONTAINER_BASE_HPP__GFLOW__