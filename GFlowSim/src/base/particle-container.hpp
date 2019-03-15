#ifndef __PARTICLE_CONTAINER_HPP__GFLOW__
#define __PARTICLE_CONTAINER_HPP__GFLOW__

#include "../gflow.hpp"
#include <set> // For storing sets of holes in the arrays
#include <unordered_map> 

namespace GFlowSimulation {

  template<typename T, typename S> T byte_cast(S x) {
    return *reinterpret_cast<T*>(&x);
  }

  class ParticleContainer : public Base {
  public:
    //! \brief Default constructor.
    ParticleContainer(GFlow*);

    //! \brief Destructor. Cleans up the array.
    ~ParticleContainer();

    //! \brief Initialize the object. 
    //! 
    //! This should be called once all particle data fields have been added.
    virtual void initialize() override;

    //! \brief Remove all halo and ghost particles.
    virtual void post_integrate() override;

    //! \brief Reserve space for particles, extending the lengths of all arrays to the requested size.
    void reserve(int);

    //! \brief Add a default particle, the properties of this particle should be set from the outside after this is called.
    void addParticle();

    //! \brief Add a particle to the container.
    void addParticle(const RealType*, const RealType*, const RealType, const RealType, const int);

    // \brief Mark a particle for removal.
    void markForRemoval(const int);

    //! \brief Remove all the particles that need to be removed, consolidate data.
    void doParticleRemoval();

    //! \brief Exchange particles with neighboring domains.
    void exchangeParticles();

    //! \brief Do a quick sort based on the particle's positions.
    void sortParticles();

    //! \brief Update the primary particle that halo particles correspond to.
    void updateHaloParticles();

    // --- Data access

    RealType *X(int);
    RealType &X(int, int);
    RealType *V(int);
    RealType &V(int, int);
    RealType *F(int);
    RealType &F(int, int);
    RealType &Sg(int);
    RealType &Im(int);

    const RealType *X(int) const;
    const RealType *V(int) const;
    const RealType *F(int) const;
    RealType Sg(int) const;
    RealType Im(int) const;
    
    //! \brief Get the type of a particle. Since this information is encoded as a float, we cannot just return a reference to it.
    int Type(int);
    void setType(int, int);
    int Id(int);

    //! \brief Get a pointer to the start of the i-th particle's data.
    RealType *particle(int);

    //! \brief Get the size of the part of the container that may contain valid particles.
    int size();

    //! \brief Get the number of particles in the container.
    int number();

    // --- Create data entries

    void addVectorData(const string&);
    void addScalarData(const string&);
    void addIntegerData(const string&);

    // --- Get data entries

    int requestVectorData(const string&);
    int requestScalarData(const string&);
    int requestIntegerData(const string&);

    // --- Clearing functions

    //! \brief Clear all the particle's velocity vectors.
    void clearV();

    //! \brief Clear all the particle's force vectors.
    //!
    //! If there is no force data, this function does nothing.
    void clearF();

  private:

    //! \brief Resize an array, increasing its capacity by a specified amount.
    //!
    //! \param additional_space The number of additional particles the resized array should be able to hold.
    inline void resize(int);

    //! \brief Swap two particles, based on their local ids.
    inline void swap_particle(int, int);

    //! \brief The data stored by the container.
    RealType* particle_data = nullptr;

    //! \brief The number of (consectutive) entries in the particle_data per particle.
    int data_width = -1;

    //! \brief The next global id that should be assigned to a particle.
    int _next_global_id = 0;

    //! \brief Whether the container has been initialized to its current configuration.
    //!
    //! This can be undone by e.g. adding another data entry.
    bool _is_initialized = false;

    //! \brief This flag is true when objects that depend on the local ids of particles in 
    //! this container should be remade. This occurs when particles are sorted or removal is
    //! carried out.
    bool needs_remake = false;

    //! \brief The offsets for particle vectorial data.
    std::map<string, int> vector_data_offsets;
    //! \brief The offsets for particle scalar data.
    std::map<string, int> scalar_data_offsets;
    //! \brief The offsets for particle integer data.
    std::map<string, int> integer_data_offsets;

    //! \brief A map between global and local particle ids, <global id, local id>.
    std::unordered_map<int, int> id_map;

    //! \brief Record where "holes" are in the particle array
    std::set<int> remove_list;

    // --- Necessary data

    //! \brief The offset for the position.
    int _pos_offset = -1;
    //! \brief The offset for the velocity.
    int _vel_offset = -1;
    //! \brief The offset for the type.
    int _type_offset = -1;
    //! \brief The offset for the global ID.
    int _id_offset = -1;

    // --- Unnecessary, but common, data

    //! \brief The offset for force.
    int _force_offset = -1;
    //! \brief The offset for the cutoff radius data.
    int _sg_offset = -1;
    //! \brief The offset for the inverse mass data.
    int _im_offset = -1;

    // --- Numbers

    //! \brief The number of particles in the container.
    //!
    //! Number is measured in particles, not float entries.
    //! This can be different than size, since there can be invalid particles (of type -1) interspersed amongst the good particles.
    int _number = 0;

    //! \brief The last part of the array that might contain valid particles.
    //!
    //! Size is measured in particles, not float entries
    //! Often, this might be the entry after the last valid particle on the processor. However, if the last valid particle was deleted, 
    //! this might not be the case. 
    int _size = 0;

    //! \brief The total capacity of the particle data arrays.
    //!
    //! Capacity is measured in particles, not float entries.
    int _capacity = 0;

    //! \brief Copy this data from force master.
    int _ntypes = -1;

  };

}
#endif // __PARTICLE_CONTAINER_HPP__GFLOW__