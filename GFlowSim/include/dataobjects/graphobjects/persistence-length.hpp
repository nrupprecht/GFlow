#ifndef __PERSISTENCE_LENGTH_HPP__GFLOW__
#define __PERSISTENCE_LENGTH_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"
#include "../../other/group.hpp"

namespace GFlowSimulation {

  class PersistenceLength : public GraphObject, public Group {
  public:
    //! \brief Default constructor
    PersistenceLength(GFlow*);

    //! \brief Constructor that sets the group.
    PersistenceLength(GFlow*, Group&);

    //! \brief Set up binning.
    void pre_integrate() override;

    //! \brief Record data.
    void post_step() override;

    //! \brief Look at the data once the simulation ends.
    void post_integrate() override;
  };

}
#endif // __PERSISTENCE_LENGTH_HPP__GFLOW__
