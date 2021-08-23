#ifndef __DOMAIN_TEST_HPP__GFLOW__
#define __DOMAIN_TEST_HPP__GFLOW__

#include "../base/domainbase.hpp"
#include "cell.hpp"

namespace GFlowSimulation {

class Domain : public DomainBase {
 public:
  //! \brief Default constructor.
  explicit Domain(GFlow *);

  //! \brief Get all the particles within a radius of another particle
  //!
  //! Fills a passed in vector with the ids of all the particles that lie within
  //! a specified distance of a given particle.
  //! This function must be overloaded by all children of DomainBase.
  void getAllWithin(int id1, std::vector<int> &neighbors, RealType distance) override;

  //! \brief Get all the particles withing a radius of some position.
  void getAllWithin(Vec X, std::vector<int> &neighbors, RealType distance) override;

  //! \brief Construct interactions for a single particle.
  //!
  //! If the insert flag is set to true, the particle will first be added to the relevent domain or data structure
  //! (if applicable), if not, it will not be added.
  //! This function should be called after construct.
  void constructFor(int id1, bool insert) override;

  //! \brief Traverse the cells structure, calling a function on pairs of particles that are close to one another.
  //!
  //! The function that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
  //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
  void traversePairs(PairFunction body) override;

  //! \brief Function that traverses ghost particles, calling a function on all ghost-real pairs of particls that are within
  //! cutoff + skin depth distance of one another.
  //!
  //! The function (a PairFunction) that is passed in should expect to receive particles' id1, id2, wrapping type (0 - no wrapping
  //! required, 1 - wrapping required), radius of particle 1, radius of particle 2, and distance between particles.
  void traverseGhostPairs(PairFunction body) override;

 private:

  // --- Overloaded functions

  //! \brief Update the linked cells structure.
  void structure_updates() override;

  //! \brief Calculates the domain cell dimensions, widths, and inverse widths given
  //! that the cutoff has been calculated.
  void calculate_domain_cell_dimensions() override;

  //! \brief Create the cells.
  void create_cells() override;

  // --- Helper functions

  //! \brief Correct a linear index for wrapping. Returns true if the index is a valid cell.
  //!
  //! If the flag is set to false, we do not wrap positions.
  inline bool correct_index(vector<int> &, bool= true);

  //! \brief Add a particle to the cell it belongs in.
  inline void add_to_cell(const RealType *, int);

  //! \brief Returns whether or not a cell is a halo cell.
  inline bool is_halo_cell(const vector<int> &);

  //! \brief A vector holding all the cells in the domain.
  std::vector<Cell> cells;

  //! \brief Function used for setting an excluded particle when calling getAllWithin
  int _exclude = -1;
};

}
#endif // __DOMAIN_TEST_HPP__GFLOW__
