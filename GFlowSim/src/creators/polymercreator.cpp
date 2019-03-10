#include "polymercreator.hpp"
// Other files
#include "../allbonded.hpp"

namespace GFlowSimulation {

  void PolymerCreator::createArea(HeadNode *head, GFlow *gflow, std::map<string, string>& variables) {
    // Get number od dimensions
    int sim_dimensions = gflow->sim_dimensions;

   

    ParseHelper parser(head);
    parser.set_variables(variables);
    parser.addValidSubheading("Number");
    parser.sortOptions();
    HeadNode *hd = nullptr;

    int number = 50;
    RealType prob = 0.2;
    RealType sg_big = 0.05;
    RealType sg_small = 0.01;
    RealType im_big = 1./sphere_volume(sg_big, sim_dimensions);
    RealType im_small = 0.05/sphere_volume(sg_small, sim_dimensions);
    RealType dx_types[2];
    int interacting_type = 0, non_interacting_type = 1;
    dx_types[interacting_type] = 0.05;
    dx_types[non_interacting_type] = 0.01;

    parser.getHeading_Optional("Number");
    hd = parser.first();
    if (hd) {
      int num;
      parser.extract_first_parameter(num, hd);
      if (num>0) number = num;
    }
    

    // The bonds object
    HarmonicBond *harmonicbonds = new HarmonicBond(gflow);
    SimData *sd = gflow->simData;
    Bounds bnds = gflow->bounds;

    // Seed global random generator
    seedNormalDistribution();

    // Keep track of the random walk
    RealType *X = new RealType[sim_dimensions], *dX = new RealType[sim_dimensions], *ZERO = new RealType[sim_dimensions];
    RealType *V = new RealType[sim_dimensions], *dV = new RealType[sim_dimensions];
    zeroVec(ZERO, sim_dimensions); 
    randomNormalVec(V, sim_dimensions); // Randomly initial orientation

    // Initialize X in the center of the bounds
    RealType dt = 0.1, dx;
    int type = 1, gid1 = -1, gid2 = -1;
    for (int d=0; d<sim_dimensions; ++d)
      X[d] = 0.5*(bnds.max[d] + bnds.min[d]);

    // Create chain
    int last_type = -1, next_type = -1;
    int counts = 0; // Links since last large ball
    for (int i=0; i<number; ++i) {
      // Swap global ids
      gid1 = gid2;
      gid2 = sd->getNextGlobalID();

      // What type of particle should be generated
      if (i==0 || drand48()<prob) next_type = interacting_type;
      else next_type = non_interacting_type;

      // Calculate dx
      if (last_type!=-1)
        dx = dx_types[next_type] + dx_types[last_type];
      else dx = 0;

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
      if (next_type==interacting_type) {
        sd->addParticle(X, ZERO, sg_big, im_big, interacting_type);
        counts = 0; // Reset counts
      }
      else {
        sd->addParticle(X, ZERO, sg_small, im_small, non_interacting_type);
        ++counts; // Increment counts
      }

      // Add bond
      if (gid1!=-1) harmonicbonds->addBond(gid1, gid2);

      // Set last type
      last_type = next_type;
    }

    // Add bonds object to gflow
    gflow->addModifier(harmonicbonds);

    // Need to relax
    Creator::hs_relax(gflow, 1.);

    RealType vsgma = 0.05;
    for (int i=0; i<sd->number(); ++i) {
      if (sd->Type(i)==interacting_type) {
        randomNormalVec(dV, sim_dimensions);
        RealType s = vsgma*randNormal();
        scalarMultVec(s, dV, sim_dimensions);
        copyVec(dV, sd->V(i), sim_dimensions);
      }
      else {
        copyVec(ZERO, sd->V(i), sim_dimensions);
      }
    }

    // Clean up
    delete [] X;
    delete [] dX;
    delete [] ZERO;
    delete [] V;
    delete [] dV;
  }

}