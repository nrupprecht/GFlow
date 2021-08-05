#ifndef __MODIFIER_HPP__
#define __MODIFIER_HPP__

#include "../gflow.hpp"
#include "simdata.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

class Modifier : public Base {
 public:
  //! \brief Default constructor.
  Modifier(GFlow *);

  //! \brief Get the remove flag.
  bool getRemove();
  //! \brief Set the remove flag.
  void setRemove(bool);

  // GFlow is a friend class
  friend class GFlow;

 protected:
  //! \brief A flag that indicates whether this modifier should be removed.
  bool remove;
};

}
#endif // __MODIFIER_HPP__
