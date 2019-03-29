#ifndef __GROUP_HPP__GFLOW__
#define __GROUP_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  /**
  *  \brief This class deals with keeping track of local and global ids of groups of particles.
  *
  *  This class DOES NOT inherit from base, so we can use multiple inheritance and not worry about diamond inheritance.
  */
  class Group {
  public:
    //! \brief Default constructor.
    Group() {};
    
    //! \brief Copy constructor
    Group(const Group&);

    //! \brief Equals
    Group& operator=(const Group &);

    // --- Accessors

    //! \brief Get the number of particles in the group.
    int size() const;

    //! \brief Get the local id of the i-th particle.
    int at(int) const;

    //! \brief Get the global id of the i-th particle.
    int g_at(int) const;

    //! \brief Find the position of the center of mass of the group.
    void findCenterOfMass(RealType*, SimData*) const;

    //! \brief Find the center of mass velocity of the group.
    void findCOMVelocity(RealType*, SimData*) const;

    //! \brief Find the net force on the group.
    void findNetForce(RealType*, SimData*) const;

    // --- Mutators

    //! \brief Add a particle, via global id, to the body.
    //!
    //! This does not check if the id is a repeat.
    virtual void add(int);

    //! \brief Update the vector of local ids to correspond to the correct particles.
    void update_local_ids(SimData*);

  protected:
    //! \brief The global ids of the particles in the group.
    vector<int> global_ids;

    //! \brief The local ids of the particles in the group.
    vector<int> local_ids;

  };

}
#endif // __GROUP_HPP__GFLOW__