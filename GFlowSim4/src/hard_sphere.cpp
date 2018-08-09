#include "hard_sphere.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Force(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION) {};

  void HardSphere::calculateForces() const {
    // Check if there are forces to calculate. Virial is reset.
    if (!Force::initCalculation()) return; 

    verletList.forceKernel(&force, &repulsion, virial);
    
    /*
    // Get the data we need
    int nverlet = verletList.vlSize(), id1(0), id2(0); // List length, id pointers
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
    int *type = Base::simData->type;

    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // Get verlet list data
    const int *verlet = verletList.getVerlet();
    RealType F[DIMENSIONS];

    // --- Go through all particles
    for (int i=0; i<nverlet; i+=2) {
      id1 = verlet[i];
      id2 = verlet[i+1];
      // Get the displacement between the particles
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
      // Check if the particles should interact
      RealType dsqr = sqr(displacement);
      if (dsqr < sqr(sg[id1] + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, displacement, normal);
        // Calculate force strength
        forceStrength(normal, distance, id1, id2);
      }
    }
    */
  }

  void HardSphere::setRepulsion(RealType r) { 
    repulsion = r; 
  }

  //! @param[in,out] normal
  //! @param[in] distance
  //! @param[in] id1
  //! @param[in] id2
  //! @param[in] simData
  //! @param[in] param_pack A parameter pack, passed in from force. Should be of the form { repulsion } (length 1).
  //! @param[in,out] virial The virial. Add the f_i \dot r_i to this.
  void HardSphere::force(RealType* normal, const RealType distance, const int id1, const int id2, const SimData *simData, 
    const RealType *param_pack, RealType &virial) 
  {
    // This should make sure that forces are zero if either object is of type -1. This does not seem to add much (any?) overhead
    RealType c1 = simData->type[id1]<0 ? 0 : 1.;
    c1 = simData->type[id2]<0 ? 0 : c1;
    // Calculate force strength
    RealType magnitude = c1*(*param_pack)*(simData->sg[id1] + simData->sg[id2] - distance);
    // Update the virial
    virial += magnitude;
    // Force strength x Normal vector -> Sets normal to be the vectorial force between the particles.
    scalarMultVec(magnitude, normal);
  }

  //! \param[in] normal The normal of the displacement between the particles. Points from id2 to id1 (?). It is safe
  //! to change this parameter.
  //! \param[in] distance The distance between the particles.
  //! \param[in] id1 The id of the first particle.
  //! \param[in] id2 The id of the second particle.
  inline void HardSphere::forceStrength(RealType *normal, const RealType distance, const int id1, const int id2) const {
    // This should make sure that forces are zero if either object is of type -1. This does not seem to add much (any?) overhead
    RealType c1 = Base::simData->type[id1]<0 ? 0 : 1.;
    c1 = Base::simData->type[id2]<0 ? 0 : c1;
    // Calculate force strength
    RealType magnitude = c1*repulsion*(simData->Sg(id1) + simData->Sg(id2) - distance);
    // Update the virial
    Force::virial += magnitude;
    // Force strength x Normal vector -> Sets normal to be the vectorial force between the particles.
    scalarMultVec(magnitude, normal);
    // Add forces to the particles.
    plusEqVec (Base::simData->f[id1], normal);
    minusEqVec(Base::simData->f[id2], normal);
  }

}
