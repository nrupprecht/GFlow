#ifndef __NEIGHBORS_HPP__
#define __NEIGHBORS_HPP__

#include "gflow.hpp"

namespace GFlowSimulation {

  class Neighbors : public Base {
    public:
      // Default constructor
      Neighbors(GFlow *);
      
      // GFlow is a friend class
      friend class GFlow;
    
    private:
      
  };

}
#endif // __NEIGHBORS_HPP__