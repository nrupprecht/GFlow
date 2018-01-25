#ifndef __INTERACTION_HPP__
#define __INTERACTION_HPP__

namespace GFlow {

  // Forward declaration
  class SimData;
  
  class Interaction {
  public:
    // Constructor
    Interaction(int a, int b) : pId1(a), pId2(b) {};

    // Calculate the force / torque from the interaction
    void doForces();

  private:
    // The particle id's of the particles in the interaction
    int pId1, pId2;

    SimData *simData;
  };

}

#endif // __INTERACTION_HPP__
