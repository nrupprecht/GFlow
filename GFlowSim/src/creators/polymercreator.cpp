#include "polymercreator.hpp"
// Other files
#include "../allbonded.hpp"
#include "../alldataobjects.hpp"
#include "../body/wallslidebody.hpp"

#include "../utility/treeparser.hpp"
#include "../modifiers/twowallmodifier.hpp"
#include "../modifiers/torque-remover.hpp"
#include "../modifiers/force-remover.hpp"

namespace GFlowSimulation {

  void PolymerCreator::createArea(HeadNode *head, GFlow *gflow, const std::map<string, string>& variables, vector<ParticleFixer>& particle_fixers) {
    // Get the dimensions
    int sim_dimensions = gflow->getSimDimensions();

    // Create parser, with variables.
    TreeParser parser(head, variables);
    // Declare valid options
    parser.addHeadingOptional("Number");
    parser.addHeadingOptional("Length");
    parser.addHeadingOptional("Phi");
    parser.addHeadingOptional("R");
    parser.addHeadingOptional("r");
    parser.addHeadingOptional("Parallel");
    parser.addHeadingOptional("IdP");
    parser.addHeadingOptional("IdC");
    parser.addHeadingOptional("DensityP");
    parser.addHeadingOptional("DensityC");
    parser.addHeadingOptional("H");
    parser.addHeadingOptional("Correlation");
    parser.addHeadingOptional("Special");
    parser.addHeadingOptional("Line");
    // Check headings for validity
    parser.check();

    // Default densities for primary and chain particles.
    RealType rhoP = 5.;
    RealType rhoC = 5.;
    imP = 1./(rhoP*sphere_volume(rP, sim_dimensions));
    imC = 1./(rhoC*sphere_volume(rC, sim_dimensions));

    // Parameters
    int number = 1;
    RealType length = 5.;
    RealType phi = 0.2;
    bool pair = false;
    RealType h = 2.5;
    useCorr = false;
    string dp = "", dc = "";
    useAngle = true;
    bool line = false;

    // Gather parameters
    parser.firstArg("Number", number);
    parser.firstArg("Length", length);
    parser.firstArg("R", rP);
    parser.firstArg("r", rC);
    parser.firstArg("Parallel", pair);
    parser.firstArg("IdP", idP);
    parser.firstArg("IdC", idC);
    parser.firstArg("DensityP", dp);
    parser.firstArg("DensityC", dc);
    parser.firstArg("Correlation", useCorr);
    parser.firstArg("Phi", phi);
    parser.firstArg("H", h);
    parser.firstArg("Line", line);

    // Potentially change masses.
    if (dp!="") {
      if (dp=="inf") imP = 0.00001; // \todo Setting this to 0 causes problems
      else imP = 1./(convert<RealType>(dp)*sphere_volume(rP, sim_dimensions));
    }
    if (dc!="") {
      if (dc=="inf") imC = 0.00001; // \todo Setting this to 0 causes problems
      else imC = 1./(convert<RealType>(dc)*sphere_volume(rC, sim_dimensions));
    }
    
    // --- Done gathering parameters, ready to act.

    // Add bonds object to gflow
    make_bond_objects(gflow);
    
    // Seed global random generator
    seedNormalDistribution();

    // Create a group correlation object
    if (correlation==nullptr && useCorr) {
      correlation = new GroupCorrelation(gflow);
      // Create a group correlation object
      correlation = new GroupCorrelation(gflow);
      correlation->setRadius(2.5*rP);
      correlation->setNBins(250);
      // Add the correlation object
      gflow->addDataObject(correlation);
    }

    if (number==2) {
      polycorr = new TwoPolymerBinForce(gflow);
      polycorr->setCType(idC);
      polycorr->setRadius(rP);
      polycorr->setMaxDistance(10*rP);
      polycorr->setNBins(250);
      gflow->addDataObject(polycorr);
    }

    if (pair) {
      createParallelPolymers(gflow, h, phi, length);
    }
    else {
      // Create all the polymers
      for (int i=0; i<number; ++i) {
        // Create the polymer, recording the group for future use.
        Group group(gflow->getSimData());
        if (line) {
          Vec Center(sim_dimensions);
          gflow->getBounds().center(Center.data);
          group = createRandomLine(gflow, Center.data, phi, length, 1.0);
        }
        else group = createRandomPolymer(gflow, length, phi, idP, idC);

        // Add a correlation item.
        if (i==0) gflow->addDataObject(new PersistenceLength(gflow, group));

        if (polycorr) {
          if (n_polymers==1) {
            polycorr->setFirstPolymer(group);
            polycorr->setFirstChain(harmonicchain);
          }
          if (n_polymers==2) {
            polycorr->setSecondPolymer(group);
            polycorr->setSecondChain(harmonicchain);
          }
        }
      }
    }
    
    // Add local fixers to the particle fixers master list
    particle_fixers.insert(particle_fixers.end(), p_fixers.begin(), p_fixers.end());
  }

  inline void PolymerCreator::make_bond_objects(GFlow *gflow) {
    // Get the dimensions
    int sim_dimensions = gflow->getSimDimensions();

    // Add bonds object to gflow
    if (sim_dimensions==2) {
      if (harmonicbonds==nullptr && !useAngle) {
        harmonicbonds = new HarmonicBond_2d(gflow);
        gflow->addBonded(harmonicbonds);
      }
      if (useAngle) harmonicchain = new AngleHarmonicChain_2d(gflow);
    }
    else if (sim_dimensions==3) {
      if (harmonicbonds==nullptr && !useAngle) {
        harmonicbonds = new HarmonicBond_3d(gflow);
        gflow->addBonded(harmonicbonds);
      }
      if (useAngle) harmonicchain = new AngleHarmonicChain_3d(gflow);
    }
    else {
      if (harmonicbonds==nullptr && !useAngle) {
        harmonicbonds = new HarmonicBond(gflow);
        gflow->addBonded(harmonicbonds);
      }
      if (useAngle) harmonicchain = nullptr; // <-----------
    }
    
    // Add the harmonic bonds modifier.
    if (harmonicbonds) {
      harmonicbonds->setSpringConstant(8.*pow(rC/0.01, sim_dimensions-1)*DEFAULT_SPRING_CONSTANT);
    }
    if (harmonicchain) {
      harmonicchain->setSpringConstant(8.*pow(rC/0.01, sim_dimensions-1)*DEFAULT_SPRING_CONSTANT);
      harmonicchain->setAngleConstant(0.05);
      gflow->addBonded(harmonicchain);
    }
  }

  void PolymerCreator::createPolymerArrangement(vector<bool>& chain_ordering, RealType phi, RealType length) {
    // Clear out vector
    chain_ordering.clear();

    // Calculate numbers of balls, spacings
    int np = ceil(0.5*phi*length/rP);

    int extra = uniform_spacing ? static_cast<int>(floor((rP-rC)/rC)) : 0;
    spacing_factor = rP / (static_cast<RealType>(extra+1)*rC);

    // Create marks
    vector<RealType> marks;
    marks.push_back(0);
    for (int i=0; i<np; ++i) marks.push_back(drand48()*length*(1.-phi));
    marks.push_back(length*(1.-phi));
    std::sort(marks.begin(), marks.end());

    // Create ordering of particle type for the polymer
    int i;
    for (i=1; i<marks.size()-1; ++i) {
      // The number of chain links to add.
      RealType expected = 0.5*(marks[i] - marks[i-1])/rC;
      // If using uniform spacing, there will (could) be chain particles "inside" the primary particles
      if (uniform_spacing) expected += extra;
      int ncl = floor(expected);

      if (drand48()<(static_cast<RealType>(ncl)-expected)) ++ncl;
      // Don't have a dangling initial chain.
      if (i==1) ncl = max(0, ncl-1);

      // Add the chain links
      for (int j=0; j<ncl; ++j) chain_ordering.push_back(false);
      // Add the primary ball
      chain_ordering.push_back(true);
    }
    // The last links
    RealType expected = 0.5*(marks[i] - marks[i-1])/rC;
    int ncl = floor(expected);
    if (drand48()<(static_cast<RealType>(ncl)-expected)) ++ncl;
    // Add the chain links

    for (int j=0; j<ncl; ++j) chain_ordering.push_back(false);
    // We have created the chain ordering, and are done.
  }

  Group PolymerCreator::createRandomPolymer(GFlow *gflow, RealType length, RealType phi, int idP, int idC) {
    // Valid probability, number, radii, and types are needed
    if (phi>1 || phi<0 || length<=0 || rP<0 || rC<0 || idP<0 || idP>=gflow->getNTypes() || idC<0 || idC>=gflow->getNTypes()) return Group(gflow);

    // Get number of dimensions
    int sim_dimensions = gflow->sim_dimensions;

    // Create a random polymer according to the specification
    RealType v_sigma = 0.1;

    // Create bonded force objects
    make_bond_objects(gflow);

    // Create a random polymer arrangement
    vector<bool> chain_ordering;
    createPolymerArrangement(chain_ordering, phi, length);

    // Random starting point
    Vec X(sim_dimensions);
    gflow->bounds.randomPoint(X.data);

    // Random initial direction
    Vec V(sim_dimensions);
    randomNormalVec(V.data, sim_dimensions);

    // Set parameters for creating the polymer
    //*****
    RealType middleX = 0.5*gflow->getBounds().wd(0)+gflow->getBounds().min[0];
    RealType dir = n_polymers & 0x1 ? 1. : -1.;
    X[0] = middleX + 0.25*dir*floor(n_polymers+1)/2;
    X[1] = 0.5*gflow->getBounds().wd(1)+gflow->getBounds().min[1] - 0.5*length;
    V[0] = 0.; V[1] = 1.;
    v_sigma = 0;
    //*****

    // Create a polymer from the specifications
    Group group = createSinglePolymer(gflow, X.data, V.data, chain_ordering, v_sigma, idP, idC);

    // Return the group
    return group;
  }

  Group PolymerCreator::createSinglePolymer(GFlow *gflow, const RealType *x0, const RealType *v0, const vector<bool>& chain_ordering, RealType sigma_v, int idP, int idC) {
    // Get simdata and bounds
    auto sd = gflow->simData;
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
    Vec X(sim_dimensions), dX(sim_dimensions), ZERO(sim_dimensions), V(sim_dimensions), dV(sim_dimensions);

    // Initialize X to start at a random position.
    RealType variance = 0.1, dx;
    int type = 1, gid1 = -1, gid2 = -1;

    // Starting point
    copyVec(x0, X);
    // Starting orientation
    copyVec(v0, V);

    // Create a group
    Group group(gflow);
    int start_id = sd->getNextGlobalID(), end_id;

    // Create chain
    for (int i=0; i<chain_ordering.size(); ++i) {
      // Swap global ids
      gid1 = gid2;
      gid2 = end_id = sd->getNextGlobalID();

      // Add gid2 to the correlation object
      if (correlation && useCorr) correlation->add(gid2);

      // What type of particle should be generated
      next_type = chain_ordering[i] ? 0 : 1;

      // Calculate dx
      if (last_type!=-1) {
        if (uniform_spacing) dx = 2*rC*spacing_factor;
        else dx = dx_types[next_type] + dx_types[last_type];
      }
      else dx = 0;

      // Advance random path
      randomNormalVec(dV.data, sim_dimensions);
      plusEqVecScaled(V.data, dV.data, sigma_v, sim_dimensions);     
      V.normalize();
      plusEqVecScaled(X.data, V.data, dx, sim_dimensions);

      // Wrap X
      for (int d=0; d<sim_dimensions; ++d) {
        if      (X[d]< bnds.min[d]) X[d] += bnds.wd(d);
        else if (X[d]>=bnds.max[d]) X[d] -= bnds.wd(d);
      }
      // Primary type
      if (next_type==0) sd->addParticle(X.data, ZERO.data, rP, imP, idP);
      // Chain link type
      else sd->addParticle(X.data, ZERO.data, rC, imC, idC);

      // Add bond - if we are using angle harmonic chains, then we only need to add particles to those.
      if (gid1!=-1 && harmonicbonds && harmonicchain!=nullptr) harmonicbonds->addBond(gid1, gid2);
      if (harmonicchain) harmonicchain->addAtom(gid2);

      // Primary type
      group.add(gid2);

      // Set last type
      last_type = next_type;
    }

    // Give the polymer some initial energy.
    RealType vsgma = 0.; //0.001;
    for (int id=start_id; id<=end_id; ++id) {
      // Random velocity direction
      randomNormalVec(V.data, sim_dimensions);
      // Random normal amount of kinetic energy
      RealType K = vsgma*fabs(randNormal());
      // Compute which velocity this implies - the global id should be the local id at this point ... hopefully
      RealType v = sqrt(2*sd->Im(id)*K);
      // Set the velocity vector of the particle
      scalarMultVec(v, V.data, sim_dimensions);
      // Add to fixers
      ParticleFixer fix = ParticleFixer(sim_dimensions, id);
      fix.velocity = V;
      p_fixers.push_back(fix);
    }

    // Increment polymer counter
    ++n_polymers;

    // Return the group
    return group;
  }

  Group PolymerCreator::createRandomLine(GFlow *gflow, const RealType *x, const RealType phi, const RealType length, const RealType spacing) {
    // Get simdata and bounds
    auto sd = gflow->simData;
    Bounds bnds = gflow->bounds;

    // Get number of dimensions
    int sim_dimensions = gflow->sim_dimensions;

    // Calculate the number of particles.
    int n = phi * length / (2 * spacing * rP);

    vector<RealType> marks;
    marks.push_back(0);
    for (int i=0; i<n; ++i) marks.push_back(drand48()*length*(1.-phi));
    std::sort(marks.begin(), marks.end());

    // Create a group
    Group group(gflow);

    // Vectors
    Vec X(sim_dimensions), ZERO(sim_dimensions);
    copyVec(x, X);
    // Add particles
    for (int i=0; i<n; ++i) {
      // Next particle position
      if (i==0)
        X[1] += (marks[i+1] - marks[i] + spacing*rP);
      else 
        X[1] += (marks[i+1] - marks[i] + 2*spacing*rP);
      // Wrap X
      for (int d=0; d<sim_dimensions; ++d) {
        if      (X[d]< bnds.min[d]) X[d] += bnds.wd(d);
        else if (X[d]>=bnds.max[d]) X[d] -= bnds.wd(d);
      }
      // Add particle
      int gid = sd->getNextGlobalID();
      sd->addParticle(X.data, ZERO.data, rP, imP, idP);

      // Add gid2 to the correlation object
      if (correlation && useCorr) correlation->add(gid);

      // Add particle
      group.add(gid);
    }

    // Increment polymer counter
    ++n_polymers;

    // Return group
    return group;
  }

  Group PolymerCreator::createOrderedLine(GFlow *gflow, const RealType *x, const RealType phi, const RealType length, const RealType spacing) {
    // Get simdata and bounds
    auto sd = gflow->simData;
    Bounds bnds = gflow->bounds;

    // Get number of dimensions
    int sim_dimensions = gflow->sim_dimensions;

    // Calculate the number of particles.
    int n = phi * length / (2 * spacing * rP);
    // Calculate distance between adjacent spheres.
    RealType distance = length / n;

    // Create a group
    Group group(gflow);

    // Vectors
    Vec X(sim_dimensions), ZERO(sim_dimensions);
    copyVec(x, X);
    // Add particles
    for (int i=0; i<n; ++i) {
      // Next particle position
      if (i==0)
        X[1] += 0.5*distance;
      else 
        X[1] += distance;
      // Wrap X
      for (int d=0; d<sim_dimensions; ++d) {
        if      (X[d]< bnds.min[d]) X[d] += bnds.wd(d);
        else if (X[d]>=bnds.max[d]) X[d] -= bnds.wd(d);
      }
      // Add particle
      int gid = sd->getNextGlobalID();
      sd->addParticle(X.data, ZERO.data, rP, imP, idP);

      // Add gid2 to the correlation object
      if (correlation && useCorr) correlation->add(gid);

      // Add particle
      group.add(gid);
    }

    // Increment polymer counter
    ++n_polymers;

    // Return group
    return group;
  }

  void PolymerCreator::createParallelPolymers(GFlow *gflow, const RealType h, const RealType phi, const RealType length) {
    // Get number of dimensions
    int sim_dimensions = gflow->sim_dimensions;

    RealType spacing = 1.;
    // Change the inverse mass to preserve uniform mass per unit length.
    imP /= spacing;

    // Create a random polymer according to the specification.
    RealType rhoP = 1.;
    RealType sigma_v = 0.;
    RealType dx = (1.+h/2)*rP; // h \el [0, 2]

    // Initial point and normal vector
    Vec x(sim_dimensions);
    // Find the center of the bounds
    gflow->bounds.center(x.data);
    x[1] -= 0.5*length;

    // Remove harmonic bonds, so the particles will not be added.
    HarmonicBond *bonds = harmonicbonds;
    harmonicbonds = nullptr;

    // Shift x to the left 
    x[0] -= dx;
    // Create first chain.
    Group group1 = createRandomLine(gflow, x.data, phi, length, spacing);
    // Shift x to the right
    x[0] += 2*dx;
    // Randomly change x[1]
    x[1] += 2*drand48()*rP;
    // Create the second chain.
    Group group2 = createRandomLine(gflow, x.data, phi, length, spacing);

    // Add the two wall modifier
    TwoWallModifier *walls = new TwoWallModifier(gflow, group1, group2);
    walls->setMinDistance(0.); // No minimum distance
    walls->setMaxDistance(5.*rP);
    walls->setBinsDataObject(100);
    walls->setMinDistanceDataObject(2.0*rP);
    walls->setMaxDistanceDataObject(4.5*rP);
    gflow->addModifier(walls);

    // Give back harmonic bonds
    harmonicbonds = bonds;
    // Correct the inverse mass.
    imP *= spacing;
  }

}