#ifndef __DEATH_HPP__GFLOW__
#define __DEATH_HPP__GFLOW__

#include "../base/modifier.hpp"

namespace GFlowSimulation {

  class DeathRate : public Modifier {
  public:
    DeathRate(GFlow*);

    DeathRate(GFlow*, const vector<RealType>&);
    
    void setRates(const vector<RealType>&);

    virtual void pre_forces() override;

  private:
    vector<RealType> deathRates;
  };

}
#endif // __DEATH_HPP__GFLOW__