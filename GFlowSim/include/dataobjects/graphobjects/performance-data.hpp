#ifndef __PERFORMANCE_DATA_HPP__GFLOW__
#define __PERFORMANCE_DATA_HPP__GFLOW__

#include "../dataobjecttypes/graphobject.hpp"

namespace GFlowSimulation {

    class PerformanceData : public GraphObject {
    public:
      PerformanceData(GFlow*);

      virtual void post_step() override;
    };

}
#endif // __PERFORMANCE_DATA_HPP__GFLOW__