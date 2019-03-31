#include "fillareacreator.hpp"

namespace GFlowSimulation {

  void FillAreaCreator::createArea(HeadNode *head, GFlow *gflow, const std::map<string, string>& variables, vector<ParticleFixer>& particle_fixers) {
    // Check if head is good
    if (head==nullptr) return;
    // For local particle templates
    std::map<string, ParticleTemplate> particle_templates;

    // Get the dimensions
    int sim_dimensions = gflow->getSimDimensions();

    // Create a parse helper
    ParseHelper parser(head);
    parser.set_variables(variables);
    // Declare valid options
    parser.addValidSubheading("Bounds");
    parser.addValidSubheading("Template");
    parser.addValidSubheading("Number");
    parser.addValidSubheading("Velocity");
    parser.addValidSubheading("Attraction");
    // Make sure only valid options were used
    if (!parser.checkValidSubHeads()) {
      cout << "Warning: Invalid Headings:\n";
      for (auto ih : parser.getInvalidSubHeads())
        cout << " -- Heading: " << ih << endl;
    }
    // Sort options
    parser.sortOptions();
    // Pointer for head nodes
    HeadNode *hd = nullptr;

    // --- Look for bounds
    parser.getHeading_Necessary("Bounds", "We need bounds!");
    Region *region = ParseConstructor::parse_region(parser.first(), variables, gflow);

    // --- Local Particle Template. Defines "types" of particles, e.g. radius distribution, density/mass, etc.
    parser.getHeading_Optional("Template");
    for (auto h : parser) ParseConstructor::parse_particle_template(h, variables, particle_templates, gflow);

    // --- Number. How to choose which particles to fill the space with.
    parser.getHeading_Necessary("Number", "We need number information!");
    // Create a structure for recording probabilities or numbers
    std::map<string, double> particle_template_numbers;
    bool useNumber = false, usePhi = false, singleType = false;
    int number(0); 
    double phi(0);
    // Get the first head with the correct heading
    hd = parser.first();
    // Check if this is using phi or number 
    if (parser.extract_parameter(hd, "Phi", phi)) usePhi = true;
    else {
      useNumber = true;
      parser.extract_first_parameter(number, hd);
    }
    // No body - we must have something in one of two forms:
    // Number: #
    // Number: Phi=#
    if (hd->subHeads.empty()) singleType = true;
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
    else {
      // Loop through subheads of hd - these are the templates
      ParseHelper subParser(hd);
      subParser.set_variables(variables);
      // Record template name, number
      for (auto m=subParser.begin(); m!=subParser.end(); ++m) {
        particle_template_numbers.insert(pair<string, double>(m.heading(), m.convert_param<double>()));
      }
    }

    // --- Velocity. How to choose particle velocities. We will find a better / more expressive way to do this later.
    parser.getHeading_Optional("Velocity");
    // Velocity option
    // 0 - Normal
    // 1 - Specified vector
    Vec Vs(sim_dimensions);
    int velocityOption = 0; // Normal velocities by default.
    hd = parser.first();
    if (hd) {
      string opt;
      parser.extract_first_parameter(opt, hd);
      if (opt=="Zero") velocityOption = 1; // Zero velocity
      else {
        parser.set_vector_argument(Vs.data, hd, sim_dimensions);
        velocityOption = 1;
      }
      // --- Other options go here
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
    
    // --- Fill with particles
    //RealType *X = new RealType[sim_dimensions], *V = new RealType[sim_dimensions];
    Vec X(sim_dimensions), V(sim_dimensions);
    RealType sigma(0.), im(0.);
    int type(0);
    // zeroVec(V, sim_dimensions);
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
        // Select a position for the particle (random uniform)
        region->pick_position(X.data);
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
          // Select a position for the particle (random uniform)
          region->pick_position(X.data);
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
  }

}