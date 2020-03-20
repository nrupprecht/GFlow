#ifndef __GROUP_HPP__GFLOW__
#define __GROUP_HPP__GFLOW__

#include "../gflow.hpp"

namespace GFlowSimulation {

  /**
  *  \brief This class deals with keeping track of local and global ids of groups of particles.
  *
  *  This class DOES NOT inherit from base, so we can use multiple inheritance and not worry about diamond inheritance.
  *  Currently, this class contains a weak pointer to a SimData object.
  */
  class Group {
  public:
    //! \brief Constructor that takes a pointer to simData.
    Group(shared_ptr<SimData>);

    //! \brief Constructor that takes a pointer to the gflow object.
    Group(GFlow*);
    
    //! \brief Copy constructor
    Group(const Group&);

    //! \brief Equals
    Group& operator=(const Group &);

    //! \brief Set this group to be a different group, erasing the information contained in the original group.
    //!
    //! The same thing as *this = <other group>.
    void set(Group&);

    // --- Accessors

    //! \brief Get the number of particles in the group.
    int size() const;

    //! \brief Whether the group is empty.
    bool empty() const;

    //! \brief Get the local id of the i-th particle.
    int at(int) const;

    //! \brief Get the global id of the i-th particle.
    int g_at(int) const;

    //! \brief Find the position of the center of mass of the group.
    void findCenterOfMass(RealType*) const;

    //! \brief Find the center of mass velocity of the group.
    void findCOMVelocity(RealType*) const;

    //! \brief Find the net force on the group.
    void findNetForce(RealType*) const;

    //! \brief Find the total mass of the objects in the group.
    RealType findTotalMass() const;

    //! \brief Find the displacement between the closest object in the group and a point.
    void findClosestObject(const RealType*, RealType*) const;

    //! \brief Add a velocity to all the particles.
    void addVelocity(RealType*) const;

    //! \brief Add whatever force is necessary to each particle to increase its acceleration by the given amount.
    void addAcceleration(RealType*) const;

    //! \brief Gets the index of where the local id is stored in the local_ids vector.
    int getIndex(int) const;

    //! \brief Check whether the group contains a particle, by local id.
    bool contains(int) const;

    //! \brief Get the sum of the distances between adjacent particles in the group, as if they formed a 
    //! chain, and we wanted to know its length.
    RealType getChainedLength() const;

    // --- Mutators

    //! \brief Add a particle, via global id, to the body.
    //!
    //! This does not check if the id is a repeat.
    virtual void add(int);

    //! \brief Check how far a position is from a group object.
    //!
    //! Depending on the type of group, this could be done in different ways.
    virtual RealType distance(const RealType*);

    //! \brief Update the vector of local ids to correspond to the correct particles.
    void update_local_ids() const;

    //! \brief Shifts all the global ids by a constant.
    void shift_global_ids(const int) const;

  protected:

    //! \brief The global ids of the particles in the group.
    mutable vector<int> global_ids;

    //! \brief The local ids of the particles in the group.
    mutable vector<int> local_ids;

    //! \brief Weak pointer to the simdata object the group exists in.
    //!
    //! I could probably replace this with a shared_ptr - it probably wouldn't matter. I will have to test and see which is faster (if it matters).
    std::shared_ptr<SimData> sim_data;

  };

}
#endif // __GROUP_HPP__GFLOW__
