#include "polymercreator.hpp"
// Other files
#include "../allbonded.hpp"

namespace GFlowSimulation {

  void PolymerCreator::createArea(HeadNode *head, GFlow *gflow, std::map<string, string>& variables) {

    int nsteps = 50;
    if (!head->params.empty()) nsteps = convert<int>(head->params[0]->partA);

    // The bonds object
    HarmonicBond *harmonicbonds = new HarmonicBond(gflow);
    int sim_dimensions = gflow->sim_dimensions;
    SimData *sd = gflow->simData;
    Bounds bnds = gflow->bounds;

    // Keep track of the random walk
    RealType *X = new RealType[sim_dimensions], *dX = new RealType[sim_dimensions], *ZERO = new RealType[sim_dimensions];
    RealType *V = new RealType[sim_dimensions], *dV = new RealType[sim_dimensions];
    zeroVec(ZERO, sim_dimensions); 
    randomNormalVec(V, sim_dimensions); // Randomly initial orientation

    // Initialize X in the center of the bounds
    RealType dt = 0.1, dx = 0.11, sg = 0.05, im = 1./sphere_volume(sg, sim_dimensions);
    int type = 1, gid1 = -1, gid2 = -1;
    for (int d=0; d<sim_dimensions; ++d)
      X[d] = 0.5*(bnds.max[d] + bnds.min[d]);

    // Create chain
    for (int i=0; i<nsteps; ++i) {
      // Swap global ids
      gid1 = gid2;
      gid2 = sd->getNextGlobalID();

      // Advance random path
      randomNormalVec(dV, sim_dimensions);
      plusEqVecScaled(V, dV, dt, sim_dimensions);      
      normalizeVec(V, sim_dimensions);
      plusEqVecScaled(X, V, dx, sim_dimensions);

      // Wrap X
      for (int d=0; d<sim_dimensions; ++d) {
        if (X[d]<bnds.min[d]) X[d] += bnds.wd(d);
        else if (X[d]>=bnds.max[d]) X[d] -= bnds.wd(d);
      }

      // Decide particle type
      type = 1;  // 1
      sg = 0.05; // If this is too small, the sectors are very small
      
      if (drand48()<0.2) {
        type = 0;
        sg = 0.05;
      }

      sd->addParticle(X, ZERO, sg, im, type);

      // Add bond
      if (gid1!=-1) harmonicbonds->addBond(gid1, gid2);
    }

    // Need to relax
    Creator::hs_relax(gflow);

    // Add bonds object to gflow
    gflow->addModifier(harmonicbonds);
    
    // Clean up
    delete [] X;
    delete [] dX;
    delete [] ZERO;
    delete [] V;
    delete [] dV;
  }

}