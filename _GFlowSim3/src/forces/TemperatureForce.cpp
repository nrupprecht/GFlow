#include "TemperatureForce.hpp"

namespace GFlow {

  TemperatureForce::TemperatureForce(RealType t) : temperature(t), viscosity(1.308e-3), radius(-1), tempDelay(default_update_delay), position(Zero), lastUpdate(-2.*default_update_delay) {};

  TemperatureForce::TemperatureForce(vec2 p, RealType r, RealType t) : temperature(t), viscosity(1.308e-3), radius(r), tempDelay(default_update_delay), position(p), lastUpdate(-2.*default_update_delay) {};
  
  TemperatureForce::TemperatureForce(RealType t, RealType v) : temperature(t), viscosity(v), radius(-1), tempDelay(default_update_delay), position(Zero), lastUpdate(-2.*default_update_delay) {};

  void TemperatureForce::_applyForce(SimData* simData) const {
    // Checks
    RealType time = simData->getTime();
    if (temperature<=0) return;
    if (time-lastUpdate<tempDelay) return;

    // Get force arrays
    RealType *fx = simData->getFxPtr();
    RealType *fy = simData->getFyPtr();
    RealType *sg = simData->getSgPtr();

    // Get the number of particles we need to update
    int domain_end = simData->getDomainEnd();

    // Precomputed values, assumes Kb = 1
    RealType DT1 = temperature/(6*viscosity*PI);
    RealType Df1 = sqrt(2.*DT1*(time-lastUpdate));

    // Add random forces
    if (radius<0)
      for (int i=0; i<domain_end; ++i) {
	RealType Df2 = sqrt(1./sg[i]);
	vec2 tForce = Df1*Df2*randNormal()*randV();
	// Add a random force
	fx[i] += tForce.x;
	fy[i] += tForce.y;
      }
    else {
      // We need the position data
      RealType *px = simData->getPxPtr();
      RealType *py = simData->getPyPtr();
      // Add temperature force to particles close enough
      for (int i=0; i<domain_end; ++i) {
	RealType rsqr = sqr(radius);
        if (sqr(simData->getDisplacement(px[i], py[i], position.x, position.y)) < rsqr) {
	  RealType Df2 = sqrt(1./sg[i]);
	  vec2 tForce = Df1*Df2*randNormal()*randV();
	  // Add a random force
	  fx[i] += tForce.x;
	  fy[i] += tForce.y;
	}
      }
    }

    // Set last update - [lastUpdate] is mutable
    lastUpdate = time;
  }

  string TemperatureForce::_summary() const {
    return ("Temperature force: T=" + toStr(temperature) + ", eta=" + toStr(viscosity)); 
  }


}
