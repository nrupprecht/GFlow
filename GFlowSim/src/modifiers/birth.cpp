#include "birth.hpp"
// Other files
#include "../base/integrator.hpp"
#include "growradius.hpp"

namespace GFlowSimulation {

  BirthRate::BirthRate(GFlow *gflow) : Modifier(gflow), minSigma(0.05) {};

  BirthRate::BirthRate(GFlow *gflow, const vector<RealType>& br) : Modifier(gflow), birthRates(br), minSigma(0.05) {};

  void BirthRate::setRate(const vector<RealType>& rates) {
    birthRates = rates;
  }

  void BirthRate::pre_forces() {
    RealType dt = Base::integrator->getTimeStep();
    
    int births = 0;
    for (int i=0; i<Base::simData->number; ++i) {

      RealType sigma = Base::simData->Sg()[i];
      if (sigma>=0.05) split(i, 5);

      /*
      int type = Base::simData->type[i];
      if (birthRates.size()<=type) continue;
      RealType rate = birthRates[type];
      RealType sigma = Base::simData->sg[i];
      // Probability is 1-exp(-rate*dt) ~ rate*dt
      if (type>-1 && sigma>=minSigma && drand48()<rate*dt) {
        // Increment births
        ++births;
        // Copy X to X1, X2
        split(i);
      }
      */
    }
    // Remake lists 
    if (births>0) Base::simData->setNeedsRemake(true);
  }

  inline void BirthRate::split(const int id, const RealType growTime) const {
    // Data structures
    RealType *Xhat = new RealType[sim_dimensions];
    RealType *dX = new RealType[sim_dimensions];
    // Set x unit vector
    zeroVec(Xhat, sim_dimensions); 
    Xhat[0] = 1.;

    // Get the parent particle's radius
    RealType rf = Base::simData->Sg()[id];
    RealType im = Base::simData->Im()[id];
    int type    = Base::simData->type[id];
    // Add a particle and a modifer for it
    int place = Base::simData->number;
    int id1 = Base::simData->getNextGlobalID();
    Base::simData->addParticle(Base::simData->X(id), Base::simData->V(id), 0.45*rf, im, type);
    Base::gflow->addModifier(new GrowRadius(Base::gflow, id1, 0.45*rf, rf, growTime));
    // Particles will be at X +/- dX
    scalarMultVec(0.5*rf, Xhat, dX, sim_dimensions);
    plusEqVec (Base::simData->X(id), dX, sim_dimensions);
    minusEqVec(Base::simData->X(place), dX, sim_dimensions);
    // Shrink first particle and put in place, add grow radius modifier
    Base::simData->Sg()[id] *= 0.45;
    int id2 = Base::simData->global_id[id];
    Base::gflow->addModifier(new GrowRadius(Base::gflow, id2, 0.45*rf, rf, growTime));

    // Clean up
    delete [] Xhat;
    delete [] dX;
  }

}