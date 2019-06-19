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

    //! \brief Whether the group is empty.
    bool empty() const;

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

    //! \brief Find the total mass of the objects in the group.
    RealType findTotalMass(SimData*) const;

    void findClosestObject(const RealType*, RealType*, SimData*) const;

    //! \brief Add a velocity to all the particles.
    void addVelocity(RealType*, SimData*) const;

    //! \brief Add whatever force is necessary to each particle to increase its acceleration by the given amount.
    void addAcceleration(RealType*, SimData*) const;

    //! \brief Gets the index of where the local id is stored in the local_ids vector.
    //!
    //! Only works if use_correspondence is true. If not, or the local id is not in the group, returns -1.
    int getIndex(int) const;

    // --- Mutators

    //! \brief Set this group to be a different group, erasing the information contained in the original group.
    //!
    //! The same thing as *this = <other group>.
    void set(Group&);

    //! \brief Add a particle, via global id, to the body.
    //!
    //! This does not check if the id is a repeat.
    virtual void add(int);

    //! \brief Update the vector of local ids to correspond to the correct particles.
    void update_local_ids(SimData*) const;

    void setUseCorrespondence(bool);

  protected:
    //! \brief Remakes to correspondence map if use_correspondence is true.
    void redo_correspondence() const;

    //! \brief The global ids of the particles in the group.
    vector<int> global_ids;

    //! \brief The local ids of the particles in the group.
    mutable vector<int> local_ids;

    //! \brief Map from local ids to place in the local id vector.
    mutable std::multimap<int, int> correspondence;
    //! \brief Whether to use the correspondence object.
    bool use_correspondence = false;

  };

}
#endif // __GROUP_HPP__GFLOW__