#include "fillareacreator.hpp"
// Other files
#include "../utility/treeparser.hpp"

namespace GFlowSimulation {

  void FillAreaCreator::createArea(HeadNode *head, GFlow *gflow, const std::map<string, string>& variables, vector<ParticleFixer>& particle_fixers) {
    // Check if head is good
    if (head==nullptr) return;
    // For local particle templates
    std::map<string, ParticleTemplate> particle_templates;

    // Get the dimensions
    int sim_dimensions = gflow->getSimDimensions();

    // Create parser, with variables.
    TreeParser parser(head, variables);
    // Declare valid options.
    parser.addHeadingNecessary("Bounds", "We need bounds!");
    parser.addHeadingNecessary("Number", "We need number information!");
    parser.addHeadingOptional("Template");
    parser.addHeadingOptional("Attraction");
    parser.addHeadingOptional("Excluded");
    parser.addHeadingOptional("Velocity");
    // Check headings for validity.
    parser.check();

    // --- Look for bounds
    Region *region = ParseConstructor::parse_region(parser.getNode("Bounds"), variables, gflow);

    // --- Look for excluded region
    vector<Region*> excluded_regions;
    if (parser.begin("Excluded")) {
      do {
        excluded_regions.push_back(ParseConstructor::parse_region(parser.getNode(), variables, gflow));
      } while (parser.next());
    }

    // --- Local Particle Template. Defines "types" of particles, e.g. radius distribution, density/mass, etc.
    if (parser.begin("Template")) {
      do {
        ParseConstructor::parse_particle_template(parser.getNode(), variables, particle_templates, gflow);
      } while (parser.next());
      // Return the parser to the original level.
      parser.end();
    }

    // --- Number. How to choose which particles to fill the space with.
    // Create a structure for recording probabilities or numbers
    std::map<string, double> particle_template_numbers;
    bool useNumber = false, usePhi = false, singleType = false;
    int number(0); 
    double phi(0);
    
    // There must be a "Number" heading, since it is required.
    parser.focus("Number");
    // Check if this is using phi or number 
    if (parser.argName()=="Phi") {
      parser.val(phi);
      usePhi = true;
    }
    else {
      parser.arg(number);
      useNumber = true;
    }
    // No body - we must have something in one of two forms:
    // Number: #
    // Number: Phi=#
    if (parser.body_size()==0) singleType = true;
    // Yes body - either the total number, or total phi is given, and particles are generated
    // randomly with certain probabilities. Subheads must be in one of the two forms:
    //
    // Number: Phi=# {
    //   Template: #[prob]
    //   ...
    // }
    // <or>
    // Number: {
    //   Template: #
    //   ...
    // }
    else if (parser.begin()) { // This will be true since the body size is not zero.
      do {
        particle_template_numbers.insert(pair<string, double>(parser.heading(), parser.arg_cast<double>()));
      } while (parser.next());
      // Return the parser to the original level.
      parser.end();
    }
    parser.up();

    // --- Velocity. How to choose particle velocities. We will find a better / more expressive way to do this later.
    Vec Vs(sim_dimensions);
    int velocityOption = 0; // Normal velocities by default.
    // Velocity option
    // 0 - Normal
    // 1 - Specified vector
    if (parser.focus("Velocity")) {
      if (parser.argName()=="Zero") velocityOption = 1; // Zero velocity
      else {
        Vs = parser.argVec();
        velocityOption = 1;
      }
      // --- Other options would go here
    }

    // --- Check that we have defined a good area
    if (!usePhi && number<=0 && singleType)
      throw BadStructure("If using a single type, we need a nonzero number of particles.");

    // Get the simdata
    SimData *simData = gflow->getSimData();

    // Velocity selecting function
    auto select_velocity = [&] (RealType *V, RealType *X, RealType sigma, RealType im, int type) -> void {
      RealType vsgma = 0.25;
      // Velocity based on KE
      double ke = fabs(vsgma*normal_dist(generator));
      double velocity = sqrt(2*im*ke/127.324);
      // Random normal vector
      randomNormalVec(V, sim_dimensions);
      // Set the velocity
      scalarMultVec(velocity, V, sim_dimensions);
    };

    // A lambda for checking if a particle falls within any excluded region.
    auto excluded_contains = [&] (Vec x) {
      bool contains = false;
      for (auto r : excluded_regions) contains |= r->contains(x.data);
      return contains;
    };
    // The maximum number of attempts to make
    int max_attempts = 50;
    
    // --- Fill with particles
    Vec X(sim_dimensions), V(sim_dimensions);
    RealType sigma(0.), im(0.);
    int type(0);

    // If we are filling to a specified packing fraction
    if (usePhi) {
      // Create discrete distribution
      vector<double> probabilities;
      vector<ParticleTemplate> template_vector;

      // Map particle type to probability
      for (auto &pr : particle_template_numbers) {
        auto it = particle_templates.find(pr.first);
        if (it==particle_templates.end()) 
          throw BadStructure("An undefined particle type was encountered: "+pr.first);
        template_vector.push_back(it->second);
        probabilities.push_back(pr.second);
      }

      // A discrete distribution we use to choose which particle template to use next
      std::discrete_distribution<int> choice(probabilities.begin(), probabilities.end());
      int i(0);
      RealType vol = 0, Vol = region->vol();
      while (vol/Vol < phi) {
        // Select a position for the particle (random uniform) that does not fall within an excluded region
        int attempts = 0;
        do {
          region->pick_position(X.data);
          ++attempts;
        } while (excluded_contains(X) && attempts<=max_attempts);

        // Choose a type of particle to create
        int pt = choice(generator);
        ParticleTemplate &particle_creator = particle_template_numbers.empty() ? particle_templates[0] : template_vector.at(pt);
        
        // Select other characteristics
        particle_creator.createParticle(X.data, sigma, im, type, i, sim_dimensions);
        // Get next global id
        int gid = simData->getNextGlobalID();
        // Add the particle
        simData->addParticle(X.data, V.data, sigma, im, type);
        // Pick an initial velocity and create the particle fixer
        if      (im==0)             V.zero();
        else if (velocityOption==0) select_velocity(V.data, X.data, sigma, im, type);
        else if (velocityOption==1) V = Vs;
        else zeroVec(V.data, sim_dimensions);
        ParticleFixer pfix(sim_dimensions, gid);
        copyVec(V.data, pfix.velocity, sim_dimensions);
        particle_fixers.push_back(pfix);
        // Increment volume and counter
        vol += sphere_volume(sigma, sim_dimensions);
        ++i;
      }

    }
    // Else, we are filling to a specified number
    else {
      // Insert the requested number of each particle type
      for (auto &pr : particle_template_numbers) {
        auto it = particle_templates.find(pr.first);
        if (it==particle_templates.end())
          throw BadStructure("An undefined particle type was encountered: "+pr.first);
        int num = static_cast<int>(pr.second);

        ParticleTemplate &particle_creator = it->second;
        for (int i=0; i<num; ++i) {
          // Select a position for the particle (random uniform) that does not fall within an excluded region
          int attempts = 0;
          do {
            region->pick_position(X.data);
            ++attempts;
          } while (excluded_contains(X) && attempts<=max_attempts);
          
          // Select other characteristics
          particle_creator.createParticle(X.data, sigma, im, type, i, sim_dimensions);
          // Get next global id
          int gid = simData->getNextGlobalID();
          // Add the particle
          simData->addParticle(X.data, V.data, sigma, im, type);
          // Pick an initial velocity and create the particle fixer
          if      (im==0)             V.zero();
          else if (velocityOption==0) select_velocity(V.data, X.data, sigma, im, type);
          else if (velocityOption==1) V = Vs;
          else zeroVec(V.data, sim_dimensions);
          ParticleFixer pfix(sim_dimensions, gid);
          copyVec(V.data, pfix.velocity, sim_dimensions);
          particle_fixers.push_back(pfix);
        }
      }
    }

    // Clean up
    delete region;
    for (auto &er : excluded_regions) delete er;
  }

}