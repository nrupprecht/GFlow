#ifndef __DEMON_HPP__GFLOW__
#define __DEMON_HPP__GFLOW__

#include "../base/modifier.hpp"
#include "../other/group.hpp"
#include "../utility/vec.hpp"
#include "../dataobjects/dataobjecttypes/graphobject.hpp"
#include "../interactions/demon_wall.hpp"

namespace GFlowSimulation {

  /**
  *  \brief A model of a Maxwell demon, which controls a door that only operates under certain conditions.
  *
  *  This class inherits from group, the group of particles are the particles that make up the ``door.''
  *  The partition is assumed to be a (d-1) dimensional plane at a certain x position.
  */
  class Demon : public Modifier, public Group {
  public:
    //! \brief Default constructor.
    Demon(GFlow*);

    //! \brief Makes sure that simdata has a side entry array, do initial checkpointing.
    virtual void pre_integrate() override;

    //! \brief Checkpointing should take place pre forces, so the forces do not have to be checkpointed
    virtual void pre_forces() override;

    //! \brief Set data in the parameters object.
    virtual void post_integrate() override;

    //! \brief The decision criteria of whether to open the door or not.
    //!
    //! \param nl The number of particles that pass from left to right.
    //! \param nr The number of particles that pass from right to left.
    //! \param el The energy that passes from left to right.
    //! \param er The energy that passes from right to left.
    virtual bool should_door_open(int, int, RealType, RealType);

    //! \brief Set the partition position.
    void setPartitionPosition(RealType);

    //! \brief Set the demon wall interaction.
    void setInteraction(shared_ptr<DemonWall>);

    //! \brief Set the start time for the door.
    void setStartTime(RealType);

    //! \brief Set the demon's tau.
    void setTau(RealType);

    //! \brief Set the type of demon.
    //!
    //! Options: 0 - simple, 1 - energy, 2 - number.
    void setDemon(int);

  private:

    //! \brief Count the particles and energy of transfered particles.
    inline void count_changes(int&, int&, RealType&, RealType&);

    //! \brief Write particles' side data to simdata.
    inline void assign_side();

    //! \brief Checkpoint the data.
    inline void checkpoint_data();

    //! \brief Reset the particles to be at the last checkpoint.
    //!
    //! The flag is whether or not to force the domains to remake.
    inline void revert_to_last_checkpoint(bool=true);

    //! \brief Open the door.
    inline void open_door();

    //! \brief Close the door.
    inline void close_door();

    //! \brief Add data to the vectors.
    inline void add_data(int, RealType);

    //! \brief Get what side a particle is on based only on position.
    inline int get_side(RealType*);

    //! \brief A function pointer to the door open/close choice function.
    bool (*check_choice) (int, int, RealType, RealType);
    //! \brief Integer representation of demon type.
    int demon_type = 0;

    //! \brief Choice function to make a direction demon.
    static bool direction_demon(int, int, RealType, RealType);
    //! \brief Choice function to make a energy demon.
    static bool energy_demon(int, int, RealType, RealType);
    //! \brief Choice function to make a number demon.
    static bool number_demon(int, int, RealType, RealType);

    //! \brief The demon wall interaction.
    shared_ptr<DemonWall> demon_interaction;

    //! \brief Whether the door is open.
    bool door_open = true;

    //! \brief Whether the door was opened last time step.
    bool door_just_closed = false;

    //! \brief When the demon should start working.
    RealType startTime = 0.;

    //! \brief The demon's time delay.
    RealType tau = 0.1;

    //! \brief The last door open/close point.
    RealType last_check = 0.;

    //! \brief Save all particles's positions. Save as a single (D*number)-dimensional array.
    vector<Vec> checkpoint_x;

    //! \brief Save all particle's velocities. Save as a single (D*number)-dimensional array.
    vector<Vec> checkpoint_v;

    //! \brief The entry for what side particles were on as of the last checkpointing.
    int side_entry = -1;

    //! \brief The x position of the partition plane.
    RealType partition_position = 0.;

    //! \brief Door half width
    RealType hwidth = 0.;

    //! \brief The positions of the particles that make up the door.
    vector<Vec> door_positions;

    // --- Statistics

    //! \brief Last energy on the left side.
    RealType El = -1.;
    //! \brief Last energy on the right side.
    RealType Er = -1.;
    //! \brief Last number of particles on the left side.
    int Nl = -1;
    //! \brief Last number of particles on the right side.
    int Nr = -1;

    int opened_tau = 0, closed_tau = 0;

    //! \brief Data objects for recording demon statistics.
    shared_ptr<GraphObject> kineticL, kineticR, numberL, numberR, current_E, current_N;

    //! \brief A pointer to a position data object that we try to find in datamaster.
    //!
    //! Having this allows us to modify the animation object to correctly animate what is happening.
    shared_ptr<class PositionData> animate_object;

    //! \brief A pointer to a parameters object.
    shared_ptr<class Parameters> parameters;

    //! \brief Save the last record time of the animation object.
    RealType animate_last_recording = 0;
    //! \brief Save the last size of the animation object
    int animate_last_size = 0;

  };

}

#endif // __DEMON_HPP__GFLOW__
