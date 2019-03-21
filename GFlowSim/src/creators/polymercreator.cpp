#include "polymercreator.hpp"
// Other files
#include "../allbonded.hpp"
#include "../alldataobjects.hpp"

namespace GFlowSimulation {

  void PolymerCreator::createArea(HeadNode *head, GFlow *gflow, std::map<string, string>& variables) {

    //createLine(head, gflow, variables);
    //return;

    

    ParseHelper parser(head);
    parser.set_variables(variables);
    parser.addValidSubheading("Number");
    parser.sortOptions();
    HeadNode *hd = nullptr;


    // Seed global random generator
    seedNormalDistribution();

    int number = 50;
    RealType prob = 0.2;
    RealType sg_big = 0.05;
    RealType sg_small = 0.01;

    parser.getHeading_Optional("Number");
    hd = parser.first();
    if (hd) {
      int num;
      parser.extract_first_parameter(num, hd);
      if (num>0) number = num;
    }
    
    // Create a group correlation object
    if (correlation==nullptr) correlation = new GroupCorrelation(gflow);
    correlation->setRadius(2.5*sg_big);

    // Seed global random generator
    seedNormalDistribution();

    createPolymer(gflow, number, prob, sg_big, sg_small, 0, 1);
  }

  void PolymerCreator::createPolymer(GFlow *gflow, int number, RealType prob, RealType sigmaP, RealType sigmaC, int idP, int idC) {
    // Valid probability, number, radii, and types are needed
    if (prob>1 || prob<0 || number<=0 || sigmaP<0 || sigmaC<0 || idP<0 || idP>=gflow->getNTypes() || idC<0 || idC>=gflow->getNTypes()) return;

    // Get number of dimensions
    int sim_dimensions = gflow->sim_dimensions;

    // Create a random polymer according to the specification
    RealType imP = 1./sphere_volume(sigmaP, sim_dimensions);
    RealType imC = 1./sphere_volume(sigmaC, sim_dimensions);
    RealType dx_types[2];
    dx_types[0] = sigmaP;
    dx_types[1] = sigmaC;
    
    // Create a group correlation object
    GroupCorrelation *correlation = new GroupCorrelation(gflow);
    correlation->setRadius(2.5*max(sigmaC, sigmaP));

    // The bonds object
    HarmonicBond *harmonicbonds = new HarmonicBond(gflow);
    SimData *sd = gflow->simData;
    Bounds bnds = gflow->bounds;

    // Keep track of the random walk
    RealType *X = new RealType[sim_dimensions], *dX = new RealType[sim_dimensions], *ZERO = new RealType[sim_dimensions];
    RealType *V = new RealType[sim_dimensions], *dV = new RealType[sim_dimensions];
    zeroVec(ZERO, sim_dimensions); 
    randomNormalVec(V, sim_dimensions); // Randomly initial orientation

    // Initialize X to start at a random position.
    RealType dt = 0.1, dx;
    int type = 1, gid1 = -1, gid2 = -1;
    for (int d=0; d<sim_dimensions; ++d)
      X[d] = drand48()*bnds.wd(d) + bnds.min[d];

    // Create chain
    int last_type = -1, next_type = -1;
    int counts = 0; // Links since last large ball
    for (int i=0; i<number; ++i) {
      // Swap global ids
      gid1 = gid2;
      gid2 = sd->getNextGlobalID();

      // Add gid2 to the correlation object
      correlation->addToGroup(gid2);

      // What type of particle should be generated
      if (i==0 || drand48()<prob) next_type = idP; // Primary type
      else next_type = idC; // Chain link type

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
      if (next_type==idP) { // Primary type
        sd->addParticle(X, ZERO, sigmaP, imP, idP);
        counts = 0; // Reset counts
      }
      else { // Chain link type
        sd->addParticle(X, ZERO, sigmaC, imC, idC);
        ++counts; // Increment counts
      }

      // Add bond
      if (gid1!=-1) harmonicbonds->addBond(gid1, gid2);

      // Set last type
      last_type = next_type;
    }

    // Add bonds object to gflow
    gflow->addModifier(harmonicbonds);
    // Add the correlation object
    gflow->addDataObject(correlation);

    // Need to relax
    Creator::hs_relax(gflow, 0.5);

    RealType vsgma = 0.001;
    for (int i=0; i<sd->size(); ++i) {
      //if (sd->Type(i)!=idC) {
        // Random velocity direction
        randomNormalVec(dV, sim_dimensions);
        // Random normal amount of kinetic energy
        RealType K = vsgma*fabs(randNormal());
        // Compute which velocity this implies
        RealType v = sqrt(2*sd->Im(i)*K);
        // Set the velocity vector of the particle
        scalarMultVec(v, dV, sim_dimensions);
        copyVec(dV, sd->V(i), sim_dimensions);
        /*
      }
      else copyVec(ZERO, sd->V(i), sim_dimensions);
      */
    }

    // Clean up
    delete [] X;
    delete [] dX;
    delete [] ZERO;
    delete [] V;
    delete [] dV;
  }

  void PolymerCreator::createLine(HeadNode *head, GFlow *gflow, std::map<string, string>& variable) {
    // Get number of dimensions
    int sim_dimensions = gflow->sim_dimensions;
    // Get the bounds from gflow
    Bounds bnds = gflow->getBounds();
    RealType width = bnds.wd(0);
    RealType phi = 0.5;
    RealType sigma = 0.05;
    // Calculate how many particles we should have
    int number = static_cast<int>(phi*width/(2*sigma));

    // Create place and offset lists
    vector<RealType> places;
    for (int i=0; i<number+1; ++i) {
      places.push_back(drand48() * width*phi);
    }
    std::sort(places.begin(), places.end());
    vector<RealType> offsets;
    for (int i=0; i<number; ++i) {
      offsets.push_back(places[i+1] - places[i]);
    }

    // Create a line of particles.
    RealType *X = new RealType[sim_dimensions], *ZERO = new RealType[sim_dimensions];
    zeroVec(X, sim_dimensions);
    zeroVec(ZERO, sim_dimensions); 
    // X will be in the center of the bounds.
    for (int d=1; d<sim_dimensions; ++d) X[d] = bnds.min[d] + 0.5*bnds.wd(d);
    SimData *sd = gflow->simData;
    // Track the x coordinate
    RealType x = offsets[0];
    // Place all balls
    for (int i=0; i<number; ++i) {
      // Set x coordinate
      X[0] = x;
      // Add particle
      sd->addParticle(X, ZERO, sigma, 0., 1); // IM is 0, and type is 1
      // Update x
      if (i<number-1) x += (offsets[i+1] + 2*sigma);
    }

    // Create and add a center correlation object
    CenterCorrelation *correlation = new CenterCorrelation(gflow);
    correlation->setRadius(2.5*sigma);
    gflow->addDataObject(correlation);

  }

}