#ifndef __FORCE_MASTER_HPP__GFLOW__
#define __FORCE_MASTER_HPP__GFLOW__

#include "../gflow.hpp"
#include "../other/timedobject.hpp"

namespace GFlowSimulation {

/*
*  \brief Force master is in charge of short range, non-bonded interaction. It keeps a record of which
*  interaction happens between pairs of particles. This allows, for example, interactions 0<-->0 to be
*  hard spheres, 0<-->1 is sphere-triangle, etc... ForceMaster is not responsible for deleting or managing
*  force objects, GFlow is.
*
*/
class ForceMaster : public Base, public TimedObject {
 public:
  //! \brief Default constructor.
  ForceMaster(GFlow *);

  //! \brief Constructor - also includes number of types.
  ForceMaster(GFlow *, int);

  virtual void initialize() override;

  //! \brief Initialize the does-interact array.
  virtual void pre_integrate() override;

  //! \brief Run all the interactions.
  void interact();

  //! \brief Run all the interactions involving ghosts.
  void interact_ghosts();

  //! \brief Get a pointer to the force that the particle pair belongs in. Null means no force.
  shared_ptr<Interaction> getInteraction(int, int);

  //! \brief Clear all the verlet lists of all the forces.
  void clear();

  //! \brief Close the interactionhandlers for all the forces.
  void close();

  // --- Accessors

  //! \brief Get the number of types of particles in the simulation.
  int getNTypes() const;

  //! \brief Get the maximum cutoff factor for any interaction that a given particle type participates in.
  RealType getMaxCutoff(int) const;

  //! \brief Get the array of max cutoffs
  const vector<RealType> &getMaxCutoff() const;

  //! \brief Get the sum of all the potential energies of all the interactions.
  RealType getTotalPotentialEnergy() const;

  //! \brief Get the sum of all the virials of all the interactions.
  RealType getTotalVirial() const;

  //! \brief Returns whether any particle types have interactions.
  bool interactingTypes() const;

  //! \brief Returns whether the particle type interacts with other particles.
  bool typeInteracts(int) const;

  //! \brief Return whether two particle types interact.
  bool typesInteract(int, int) const;

  // --- Mutators

  //! \brief Set the number of particle types.
  void setNTypes(int);

  //! \brief Set the force in the force grid - this also adds it to the force vector here and in the GFlow
  //! object if it is not already in those locations.
  void setInteraction(int, int, const shared_ptr<Interaction> &, bool= true);

  //! \brief Set all the interactions to be of the same type.
  void setInteraction(const shared_ptr<Interaction> &);

  //! \brief Set the calculate potential flag for all interactions.
  void setCalculatePotential(bool);

  //! \brief Set the calculate virial flag for all interactions.
  void setCalculateVirial(bool);

  //! \brief Initialize the does interact array.
  void initialize_does_interact();

  // Interaction handler is a friend class.
  friend class InteractionHandler;

 private:
  //! \brief Ask all the interactions what their timescales are.
  inline void compute_timescale();

  //! \brief A timescale for the forces. The max timestep should be some fraction of this number.
  RealType time_scale = -1;

  //! \brief The fraction of the min timescale that we should use as the max timestep.
  RealType time_scale_factor = 0.02;

  //! \brief Interaction grid.
  vector<vector<shared_ptr<Interaction> > > grid;

  //! \brief Cutoff factors for each particle type
  vector<RealType> max_cutoffs;

  //! \brief Pointers to all the forces that exist in the simulation.
  vector<shared_ptr<Interaction> > interactions;

  //! \brief Number of particle types.
  int ntypes;

  //! \brief Whether the i-th particle type interacts with any particles.
  vector<bool> doesInteract;

  //! \brief Whether there are any interactions.
  bool any_interactions = false;
};

}
#endif // __FORCE_MASTER_HPP__GFLOW__