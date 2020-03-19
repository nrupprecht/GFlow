#include "fillareacreator.hpp"
// Other files
#include "../utility/treeparser.hpp"

#include <functional>

namespace GFlowSimulation {

  void FillAreaCreator::createArea(HeadNode *head, GFlow *gflow, const std::map<string, string>& variables, vector<ParticleFixer>& particle_fixers) {
    // Check if head is good
    if (head==nullptr) return;

    // Get the dimensions
    int sim_dimensions = gflow->getSimDimensions();

    // Get the number of particle fixers - any additional elements were added by this creator.
    int pf_size = particle_fixers.size();

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
    bool singleType = false;
    // Number option: 0 - phi, 1 - rho, 2 - number
    int option = 0;
    int number(0); 
    double phi(0);
    
    // There must be a "Number" heading, since it is required.
    parser.focus0("Number");
    // Volume density
    if (parser.argName()=="Phi") {
      parser.val(phi);
      option = 0;
    }
    // Number density
    else if (parser.argName()=="Rho") {
      parser.val(phi);
      option = 1;
    }
    else {
      parser.arg(number);
      option = 2;
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

    // --- Velocity. How to choose particle velocities. We will find a better / more expressive way to do this later.
    Vec Vs(sim_dimensions);
    int velocityOption = 2;         // Use temperature by default.
    RealType kinetic = 0.00196;     // Default kinetic energy factor --> 0.25/127.324
    RealType temperature = 0.00156; // Default temperature
    // Velocity option
    // 0 - Normal
    // 1 - Specified vector
    // 2 - Exponential distribution of kinetic energies.
    if (parser.focus0("Velocity")) {
      // Choose velocity to be a random normal magnitude vector
      if (parser.argName()=="" || parser.argName()=="Normal") {
        // Get the target kinetic energy, if specified.
        parser.val(kinetic);
        // Set velocity option
        velocityOption = 0;
      }
      else if (parser.argName()=="Temperature") {
        parser.val(temperature);	
        // Set velocity option
        velocityOption = 2;
      }
      // Otherwise, a specific velocity vector is specified.
      else if (parser.argName()=="Zero") velocityOption = 1; // Zero velocity
      else { 
        // Get the velocity vector
        Vs = parser.argVec();
        // Set velocity option
        velocityOption = 1;
      }
      // Initialize particles to a certain temperature

      // --- Other options would go here
    }

    // --- Check that we have defined a good area
    if (option!=0 && number<=0 && singleType)
      throw BadStructure("If using a single type, we need a nonzero number of particles.");

    // Get the simdata
    auto simData = gflow->getSimData();

    // Normal velocities
    auto choose_by_normal = [&] (RealType *V, RealType *X, RealType sigma, RealType im, int type) -> void {
      // Velocity based on KE
      double ke = fabs(normal_dist(generator));
      double velocity = sqrt(2*im*ke*kinetic);
      // Random normal vector
      randomNormalVec(V, sim_dimensions);
      // Set the velocity
      scalarMultVec(velocity, V, sim_dimensions);
    };

    // Constant velocity vector
    auto constant_velocity = [&] (RealType *V, RealType *X, RealType sigma, RealType im, int type) -> void {
      // Velocity based on KE
      copyVec(Vs.data, V, sim_dimensions);
    };

    // Expontial distribution
    std::exponential_distribution<double> exp_distribution(1./(gflow->KB * temperature));
    auto choose_by_temperature = [&] (RealType *V, RealType *X, RealType sigma, RealType im, int type) -> void {
      double ke = exp_distribution(generator);
      double velocity = sqrt(2*im*ke);
      // Random normal vector
      randomNormalVec(V, sim_dimensions);
      // Set the velocity
      scalarMultVec(velocity, V, sim_dimensions);
    };

    // Choose the function to use to assign velocities
    std::function<void(RealType*, RealType*, RealType, RealType, int)> select_velocity;
    if      (velocityOption==0) select_velocity = choose_by_normal;
    else if (velocityOption==1) select_velocity = constant_velocity;
    else if (velocityOption==2) select_velocity = choose_by_temperature;

    // A lambda for checking if a particle falls within any excluded region.
    auto excluded_contains = [&] (Vec x) {
      bool contains = false;
      for (auto r : excluded_regions) contains |= r->contains(x.data);
      return contains;
    };
    // The maximum number of attempts to make
    int max_attempts = 50;

    
    
    // --- Fill the region with particles
    Vec X(sim_dimensions), V(sim_dimensions);
    RealType sigma(0.), im(0.);
    int type(0), N(0);
    RealType ke(0.), Vol = region->vol();
    // Correct number/volume if we are doing a parallel setup, since some particles will be initialized on another processor.
    // We can only handle this case for now, since we need to intersect the region with the process bounds.
    if (dynamic_cast<RectangularRegion*>(region)) {
      //! \todo This only works if the region is rectangular.
      auto process_bounds = gflow->getTopology()->getProcessBounds();
      auto region_bounds = region->get_bounding_box();
      auto intersection = Bounds::intersection(process_bounds, region_bounds);
      // Correct number, if necessary.
      number = number * region_bounds.vol() / region->vol();
      Vol = intersection.vol();
      // Make the region into the intersection region.
      delete region;
      region = new RectangularRegion(sim_dimensions, intersection);
    }

    // If we are filling to a specified volume fraction or number density.
    if (option==0 || option==1) {
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

      // Keep track of number and total volume of particle as we add them.
      RealType vol = 0;
      // A function that checks if we should keep adding particles.
      auto keep_adding = [&] () -> bool {
        // Phi - fill to a volume fraction.
        if (option==0) return vol/Vol<phi;
        // Rho - fill to a number density.
        else if (option==1) return N/Vol<phi;
        // Else, error
        else return false;
      };

      // Add particles as long as we should
      while (keep_adding()) {
        // If a large amount of the region that a processor manages is excluded, then that processor creates particles
        // at much higher densities than it should. We do the following so that the density of particles *in the nonexcluded region* 
        // is the specified density.
        // We pick a position, and then count the volume, number, and energy of the particle whether or not we actually add the particle
        // to the simulation.
        region->pick_position(X.data);

        // Choose a type of particle to create
        int pt = choice(generator);
        ParticleTemplate &particle_creator = particle_template_numbers.empty() ? particle_templates[0] : template_vector.at(pt);
        
        // Select other characteristics
        particle_creator.createParticle(X.data, sigma, im, type, N, sim_dimensions);
        // Pick an initial velocity and create the particle fixer
        if (im==0) V.zero();
        else select_velocity(V.data, X.data, sigma, im, type);

        // Add the particle - only if not in the excluded region, see above note.
        if (!excluded_contains(X)) {
          // Get next global id
        int gid = simData->getNextGlobalID();
          simData->addParticle(X.data, V.data, sigma, im, type);
          ParticleFixer pfix(sim_dimensions, gid);
          pfix.velocity = V;
          particle_fixers.push_back(pfix);
        }
        // Increment volume, KE, and counter
        vol += sphere_volume(sigma, sim_dimensions);
        ke  += 0.5*sqr(V)/im;
        ++N;
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
        for (; N<num; ++N) {
          // Select a position for the particle (random uniform) that does not fall within an excluded region
          int attempts = 0;
          do {
            region->pick_position(X.data);
            ++attempts;
          } while (excluded_contains(X) && attempts<=max_attempts);
          
          // Select other characteristics
          particle_creator.createParticle(X.data, sigma, im, type, N, sim_dimensions);
          // Get next global id
          int gid = simData->getNextGlobalID();
          // Add the particle
          simData->addParticle(X.data, V.data, sigma, im, type);
          // Pick an initial velocity and create the particle fixer
          if (im==0) V.zero();
          else select_velocity(V.data, X.data, sigma, im, type);
          ParticleFixer pfix(sim_dimensions, gid);
          pfix.velocity = V;
          particle_fixers.push_back(pfix);
          // Increment KE
          ke  += 0.5*sqr(V)/im;
        }
      }
    }

    // If using temperature based velocity initialization, make sure the temperature is exactly right.
    if (velocityOption==2) {
      RealType Teff = ke / (N/2 * sim_dimensions * gflow->getKB()); 
      RealType frac = sqrt(temperature/Teff);
      for (int i=pf_size; i<particle_fixers.size(); ++i) particle_fixers[i].velocity *= frac;
    }

    // Clean up
    delete region;
    for (auto &er : excluded_regions) delete er;
  }

}
