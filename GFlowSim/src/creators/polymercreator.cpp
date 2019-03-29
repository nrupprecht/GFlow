#include "polymercreator.hpp"
// Other files
#include "../allbonded.hpp"
#include "../alldataobjects.hpp"

namespace GFlowSimulation {

  void PolymerCreator::createArea(HeadNode *head, GFlow *gflow, std::map<string, string>& variables) {
    // Set up parser
    ParseHelper parser(head);
    parser.set_variables(variables);
    parser.addValidSubheading("Number");
    parser.addValidSubheading("Length");
    parser.addValidSubheading("Phi");
    parser.addValidSubheading("R");
    parser.addValidSubheading("r");
    parser.addValidSubheading("Correlation");
    parser.sortOptions();
    HeadNode *hd = nullptr;

    // Parameters
    int number = 1;
    RealType length = 5.;
    RealType phi = 0.2;
    bool useAngle = false;

    // Extract the number of chains to create.
    parser.getHeading_Optional("Number");
    hd = parser.first();
    if (hd) {
      int num;
      parser.extract_first_parameter(num, hd);
      if (num>0) number = num;
    }

    // Extract the number of atoms per chain.
    parser.getHeading_Optional("Length");
    hd = parser.first();
    if (hd) {
      int l;
      parser.extract_first_parameter(l, hd);
      if (l>0) length = l;
    }

    // Extract the probability that an atom is a large particle.
    parser.getHeading_Optional("Phi");
    hd = parser.first();
    if (hd) {
      RealType p;
      parser.extract_first_parameter(p, hd);
      if (p>0) phi = p;
    }

    // Extract the radius of the large particles.
    parser.getHeading_Optional("R");
    hd = parser.first();
    if (hd) {
      RealType R;
      parser.extract_first_parameter(R, hd);
      if (R>0) rP = R;
    }

    // Extract the radius of the small (chain link) particles.
    parser.getHeading_Optional("r");
    hd = parser.first();
    if (hd) {
      RealType r;
      parser.extract_first_parameter(r, hd);
      if (r>0) rC = r;
    }

    // Whether to use a correlation functino or not.
    parser.getHeading_Optional("Correlation");
    hd = parser.first();
    if (hd) {
      bool c;
      parser.extract_first_parameter(c, hd);
      useCorr = c;
    }

    // --- Done gathering parameters, ready to act.

    // Get number of dimensions
    int sim_dimensions = gflow->sim_dimensions;
    
    // Seed global random generator
    seedNormalDistribution();

    // Create a group correlation object
    if (correlation==nullptr) {
      correlation = new GroupCorrelation(gflow);
      // Create a group correlation object
      correlation = new GroupCorrelation(gflow);
      correlation->setRadius(4.0*rP);
      correlation->setNBins(250);
      // Add the correlation object
      gflow->addDataObject(correlation);
    }

    // Add bonds object to gflow
    if (harmonicbonds==nullptr) {
      if (sim_dimensions==2) {
        harmonicbonds = new HarmonicBond_2d(gflow);
        if (useAngle) harmonicchain = new AngleHarmonicChain_2d(gflow, 0.1*DEFAULT_SPRING_CONSTANT);
      }
      else if (sim_dimensions==3) {
        harmonicbonds = new HarmonicBond_3d(gflow);
        if (useAngle) harmonicchain = new AngleHarmonicChain_3d(gflow, 0.1*DEFAULT_SPRING_CONSTANT);
      }
      else {
        harmonicbonds = new HarmonicBond(gflow);
        if (useAngle) harmonicchain = nullptr; // <-----------
      }
      // Add the harmonic bonds modifier.
      gflow->addModifier(harmonicbonds);
      gflow->addModifier(harmonicchain);
      // Adjust strength of harmonic bonds
      harmonicbonds->setSpringConstant(8.*pow(rC/0.01, sim_dimensions-1)*DEFAULT_SPRING_CONSTANT);
    }

    // Seed global random generator
    seedNormalDistribution();

    // Create all the polymers
    for (int i=0; i<number; ++i) {
      createRandomPolymer(gflow, length, phi, 0, 1);
      ++n_polymers;
    }

    // Need to relax
    Creator::hs_relax(gflow, 0.5);
    
    // Give random normal velocities.
    RealType vsgma = 0.001;
    SimData *sd = gflow->simData;
    for (int i=0; i<sd->size(); ++i) {
      // Random velocity direction
      randomNormalVec(sd->V(i), sim_dimensions);
      // Random normal amount of kinetic energy
      RealType K = vsgma*fabs(randNormal());
      // Compute which velocity this implies
      RealType v = sqrt(2*sd->Im(i)*K);
      // Set the velocity vector of the particle
      scalarMultVec(v, sd->V(i), sim_dimensions);
    }
  }

  void PolymerCreator::createRandomPolymer(GFlow *gflow, RealType length, RealType phi, int idP, int idC) {
    // Valid probability, number, radii, and types are needed
    if (phi>1 || phi<0 || length<=0 || rP<0 || rC<0 || idP<0 || idP>=gflow->getNTypes() || idC<0 || idC>=gflow->getNTypes()) return;

    // Get number of dimensions
    int sim_dimensions = gflow->sim_dimensions;

    // Create a random polymer according to the specification
    RealType rhoP = 10.;
    RealType rhoC = 10.;
    imP = 1./(rhoP*sphere_volume(rP, sim_dimensions));
    imC = 1./(rhoC*sphere_volume(rC, sim_dimensions));

    // Calculate numbers of balls, spacings
    int np = floor(0.5*phi*length/rP);
    int nc = ceil (0.5*(1.-phi)*length/rC);

    // Create marks
    vector<RealType> marks;
    marks.push_back(0);
    for (int i=0; i<np; ++i) marks.push_back(drand48()*length*(1.-phi));
    marks.push_back(length*(1.-phi));
    std::sort(marks.begin(), marks.end());

    // Create ordering of particle type for the polymer
    vector<bool> chain_ordering;
    int i;
    for (i=1; i<marks.size()-1; ++i) {
      // The number of chain links to add.
      int ncl = ceil(0.5*(marks[i] - marks[i-1])/rC);
      // Add the chain links
      for (int j=0; j<ncl; ++j) chain_ordering.push_back(false);
      // Add the primary ball
      chain_ordering.push_back(true);
    }
    // The last links
    int ncl = ceil(0.5*(marks[i] - marks[i-1])/rC);
    // Add the chain links
    for (int j=0; j<ncl; ++j) chain_ordering.push_back(false);

    // Random starting point
    Bounds bnds = gflow->bounds;
    RealType *X = new RealType[sim_dimensions];
    for (int d=0; d<sim_dimensions; ++d)
      X[d] = drand48()*bnds.wd(d) + bnds.min[d];

    // Random initial direction
    RealType *V = new RealType[sim_dimensions];
    randomNormalVec(V, sim_dimensions);

    // Create a polymer from the specifications
    createSinglePolymer(gflow, X, V, chain_ordering, 0.1, idP, idC);

    // Clean up
    delete [] X;
    delete [] V;
  }

  void PolymerCreator::createSinglePolymer(GFlow *gflow, const RealType *x0, const RealType *v0, const vector<bool>& chain_ordering, RealType sigma_v, int idP, int idC) {
    // Get simdata and bounds
    SimData *sd = gflow->simData;
    Bounds bnds = gflow->bounds;

    // Get number of dimensions
    int sim_dimensions = gflow->sim_dimensions;

    // Spacings
    RealType dx_types[2];
    dx_types[0] = rP;
    dx_types[1] = rC;
    // Counting and bookkeeping
    int last_type = -1, next_type = -1;

    // Keep track of the random walk
    RealType *X = new RealType[sim_dimensions], *dX = new RealType[sim_dimensions], *ZERO = new RealType[sim_dimensions];
    RealType *V = new RealType[sim_dimensions], *dV = new RealType[sim_dimensions];
    zeroVec(dX, sim_dimensions);
    zeroVec(ZERO, sim_dimensions);
    zeroVec(dV, sim_dimensions);

    // Initialize X to start at a random position.
    RealType variance = 0.1, dx;
    int type = 1, gid1 = -1, gid2 = -1;

    // Starting point
    copyVec(x0, X, sim_dimensions);
    // Starting orientation
    copyVec(v0, V, sim_dimensions);

    // Create chain
    
    for (int i=0; i<chain_ordering.size(); ++i) {
      // Swap global ids
      gid1 = gid2;
      gid2 = sd->getNextGlobalID();

      // Add gid2 to the correlation object
      if (correlation && useCorr) correlation->addToGroup(gid2);

      // What type of particle should be generated
      next_type = chain_ordering[i] ? 0 : 1;

      // Calculate dx
      if (last_type!=-1)
        dx = dx_types[next_type] + dx_types[last_type];
      else dx = 0;

      // Advance random path
      randomNormalVec(dV, sim_dimensions);
      plusEqVecScaled(V, dV, sigma_v, sim_dimensions);      
      normalizeVec(V, sim_dimensions);
      plusEqVecScaled(X, V, dx, sim_dimensions);

      // Wrap X
      for (int d=0; d<sim_dimensions; ++d) {
        if      (X[d]< bnds.min[d]) X[d] += bnds.wd(d);
        else if (X[d]>=bnds.max[d]) X[d] -= bnds.wd(d);
      }

      // Primary type
      if (next_type==idP) sd->addParticle(X, ZERO, rP, imP, idP);
      // Chain link type
      else sd->addParticle(X, ZERO, rC, imC, idC);

      // Add bond
      if (gid1!=-1) harmonicbonds->addBond(gid1, gid2);
      if (harmonicchain) harmonicchain->addAtom(gid2);

      // Set last type
      last_type = next_type;
    }

  }

}