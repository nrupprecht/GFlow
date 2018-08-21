#ifndef __CONSUMPTION_HPP__GFLOW__
#define __CONSUMPTION_HPP__GFLOW__

#include "../base/interaction.hpp"

namespace GFlowSimulation {

  class Consumption : public Interaction {
    //! @brief Constructor
    Consumption(GFlow *);

    //! @param[in] normal
    //! @param[in] distance
    //! @param[in] id1
    //! @param[in] id2
    //! @param[in] simData
    //! @param[in] param_pack A parameter pack, passed in from force. Contains characteristic 
    //! constants of the force, and extra data the force needs.
    //! @param[in,out] data_pack Data to be updated by the function.
    static void consume(RealType*, const RealType, const int, const int, SimData*, const RealType*, RealType*);

  };

}
#endif // __CONSUMPTION_HPP__GFLOW__