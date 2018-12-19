#include "fileparsecreator.hpp"
// Other files
#include "../utility/printingutility.hpp"
#include "../allmodifiers.hpp"

namespace GFlowSimulation {

  FileParseCreator::FileParseCreator(int argc, char **argv) : Creator(argc, argv), configFile(""), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    // Seed generators here
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  FileParseCreator::FileParseCreator(ArgParse *p) : Creator(p), configFile(""), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    // Seed generators here
    seedGenerator(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  FileParseCreator::FileParseCreator(ArgParse *p, string f) : Creator(p), configFile(f), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    // Seed generators here
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  GFlow* FileParseCreator::createSimulation() {
    // Timing 
    auto start_time = current_time();

    // Command line arguments
    RealType skinDepth = -1.;
    RealType cell_size = -1;

    if (parserPtr) {
      parserPtr->get("skinDepth", skinDepth);
      parserPtr->get("cell_size", cell_size);
    }

    // We treat the file as one giant body. Get that body.
    build_message = "Starting file parse... ";
    FileParse parser;
    HeadNode *root = nullptr;
    try {
      root = parser.parseFile(configFile); // A file parse class does the parsing
    }
    catch (FileParse::UnexpectedToken ut) {
      cout << "Caught unexpected token error from file parsing. Trace:\n";
      cout << parser.getMessage();
      throw;
    }
    build_message += "Done.\n";

    // Get the message from the parser.
    parse_message = parser.getMessage();

    // Save the file we just loaded
    string file = copyFile();

    // Sort and collect options
    build_message += "Collecting options... ";
    std::multimap<string, HeadNode*> options;
    for (const auto &op : root->subHeads) {
      options.insert(std::pair<string, HeadNode*> (op->heading, op));
    }
    build_message += "Done.\n";

    // Create the scenario from the options
    if (gflow) delete gflow;
    // Set skin depth
    //if (skinDepth>0) gflow->domain->setSkinDepth(skinDepth);
    //if (cell_size>0) gflow->domain->setCellSize(cell_size);
    // Create from the options
    build_message += "Starting simulation setup... ";
    try {
      createFromOptions(root);
    }
    catch (BadStructure bs) {
      cout << "Caught Bad Structure error: " << bs.message << endl;
      cout << "Build Message:\n" << build_message << endl;
      // Print the parse tree to a file so we can debug
      ofstream fout("ParseTrace.txt");
      if (fout.fail());
      else {
        fout << parse_message;
        fout.close();
      }
      throw;
    }
    build_message += "Done.\n";

    // Clean up
    delete root;
    
    // Timing
    auto end_time = current_time();
    gflow->dataMaster->setInitializationTime(time_span(end_time, start_time));

    gflow->dataMaster->giveFile("setup.txt", file);

    // Clean up and return
    return gflow;
  }

  inline void FileParseCreator::createFromOptions(HeadNode *head) {
    // Create a parse helper
    ParseHelper parser(head);
    // Declare valid options
    parser.addValidSubheading("Dimensions");
    parser.addValidSubheading("Seed");
    parser.addValidSubheading("Bounds");
    parser.addValidSubheading("Boundary");
    parser.addValidSubheading("Gravity");
    parser.addValidSubheading("Attraction");
    parser.addValidSubheading("Integrator");
    parser.addValidSubheading(Types_Token);
    parser.addValidSubheading(Interactions_Token);
    parser.addValidSubheading("Tempate");
    parser.addValidSubheading("Fill-area");
    parser.addValidSubheading("Particle");
    parser.addValidSubheading("Creation");
    parser.addValidSubheading("Destruction");
    parser.addValidSubheading("Modifier");
    parser.addValidSubheading("MaxDomainUpdateDelay");
    parser.addValidSubheading("MaxDT");
    parser.addValidSubheading("MinDT");
    parser.addValidSubheading("Reconcile");
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

    // --- Look for dimensions
    parser.getHeading_Optional("Dimensions");
    hd = parser.first();
    int sim_dimensions = 2;
    if (hd) {
      parser.extract_first_parameter(sim_dimensions);
      if (sim_dimensions<=0) throw BadDimension();
    }
    gflow = new GFlow(sim_dimensions);

    // --- Look for seed information for random generators
    parser.getHeading_Optional("Seed");
    hd = parser.first();
    // Read in a random seed
    if (parser.extract_first_parameter(seed));
    // Generate a random seed
    else seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    srand48(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);

    // --- Look for bounds
    parser.getHeading_Necessary("Bounds");
    hd = parser.first();
    ParseHelper subParser(hd);
    simBounds = Bounds(subParser.size());
    int d=0;
    for (int d=0; d<sim_dimensions; ++d) {
      subParser.extract_parameter(simBounds.min[d], d, 0);
      subParser.extract_parameter(simBounds.max[d], d, 1);
    }
    gflow->setBounds(simBounds);

    // --- Look for boundary conditions
    parser.getHeading_Optional("Boundary");
    hd = parser.first();
    if (hd) {
      string opt;
      if (parser.extract_first_parameter(opt)) {
        if (opt=="Repulse") gflow->setAllBCs(BCFlag::REPL);
        else if (opt=="Wrap") gflow->setAllBCs(BCFlag::WRAP);
        else if (opt=="Reflect") gflow->setAllBCs(BCFlag::REFL);
        else throw BadStructure("Unrecognized boundary condition ["+opt+"].");
      }
      else {
        ParseHelper subParser(hd);
        subParser.sortOptions();
        // Record template name, number
        int d=0;
        for (auto m=subParser.begin() ; m!=subParser.end(); ++m, ++d) {
          opt = m.convert_param<string>();
          if (opt=="Repulse") gflow->setBC(d, BCFlag::REPL);
          else if (opt=="Wrap") gflow->setBC(d, BCFlag::WRAP);
          else if (opt=="Reflect") gflow->setBC(d, BCFlag::REFL);
          else throw BadStructure("Unrecognized boundary condition ["+opt+"].");
        }
      }
    }
    else gflow->setAllBCs(BCFlag::WRAP);

    parser.getHeading_Optional("Gravity");
    hd = parser.first();
    if (hd) {
      RealType *g = new RealType[sim_dimensions];
      parser.set_vector_argument(g, hd, sim_dimensions);
      gflow->addModifier(new ConstantAcceleration(gflow, g));
      delete [] g;
    }

    // --- Attraction towards the center of the domain
    parser.getHeading_Optional("Attraction");
    RealType g;
    if (parser.extract_first_parameter(g)) gflow->setAttraction(g);

    // --- Look for integrator
    parser.getHeading_Optional("Integrator");
    hd = parser.first();
    if (hd) gflow->integrator = choose_integrator(hd);
    else gflow->integrator = new VelocityVerlet(gflow);

    // --- Look for number of particle types
    parser.getHeading_Optional(Types_Token);
    int n;
    if(parser.extract_first_parameter(NTypes));
    else NTypes = 1;
    gflow->forceMaster->setNTypes(NTypes);

    // --- Look for interactions
    parser.getHeading_Necessary(Interactions_Token); 
    hd = parser.first();   
    // If we expect something like ": Random" instead of a force grid
    if (hd->subHeads.empty()) {
      if (hd->params.size()!=1) 
        throw BadStructure("Expect one parameter for force grid, since the body is empty.");
      // Read the parameter and do what it says
      if (hd->params[0]->partA=="Random") {
        // Assign random interactions, either LennardJones or HardSphere (for now), and with equal probability (for now)
        makeRandomForces();
      }
      else throw BadStructure("Unrecognized force grid parameter.");
    }
    // Otherwise 
    else {
      // Collect all the interactions we need
      std::map<string, Interaction*> interactions;
      // The body specifies the force grid
      for (auto fg : hd->subHeads) {
        // Look for type1, type2, interaction-type [, possibly "R"]
        if (fg->params.size()!=3 && fg->params.size()!=4) 
          throw BadStructure("Force grid needs three or four parameters per line to specify interaction, found "+toStr(fg->params.size())+".");
        // Check if we don't want reflexive forces
        if (fg->params.size()==4 && fg->params[3]->partA!="NR")
          throw BadStructure("Fourth argument must be 'NR' if anything.");
        // What type of interaction do we need
        string i_token = fg->params[2]->partA;
        // Check if this is a new interaction type
        if (interactions.find(i_token)==interactions.end()) 
          interactions.insert(std::pair<string, Interaction*>(i_token, choose_interaction(fg)));
        // Get the type ids
        int t1 = convert<int>(fg->params[0]->partA), t2 = convert<int>(fg->params[1]->partA);
        // Add the interaction
        gflow->forceMaster->setInteraction(t1, t2, interactions.find(i_token)->second);
        // We already checked that the fourth parameter is an 'R' if it exists. Only add if t1!=t2
        if (fg->params.size()!=4 && t1!=t2) 
          gflow->forceMaster->setInteraction(t2, t1, interactions.find(i_token)->second);
      }
    }

    // --- Global Particle Template. Defines "types" of particles, e.g. radius distribution, density/mass, etc.
    parser.getHeading_Optional("Template");
    for (auto h : parser) getParticleTemplate(h, global_templates);

    // --- Fill an area with particles
    parser.getHeading_Optional("Fill-area");
    for (auto h : parser) fillArea(h); 

    // --- Create a single particle
    parser.getHeading_Optional("Particle");
    for (auto h : parser) createParticle(h);

    // --- Rules for creating particles
    parser.getHeading_Optional("Creation");
    hd = parser.first();
    if (hd) {
      // The subheads have no bodies, and are of the form: [type i] : [type j], [rate].
      // For now, we only allow reproduction to create the same type.
      std::map<int, RealType> birth;
      for (auto s : hd->subHeads) {
        if (s->params.size()!=2) 
          throw BadStructure("Expected a type and a rate (2 options). Found "+toStr(s->params.size())+".");
        birth.insert(std::pair<int, RealType>(
          convert<int>(s->heading), 
          convert<RealType>(s->params[1]->partA)
        ));
      }
      // Create a vector of birth rates - we must already know NTypes.
      if (!birth.empty()) {
        vector<RealType> birthRates(NTypes, 0);
        for (int i=0; i<NTypes; ++i) {
          if (contains(birth, i)) birthRates[i] = birth.find(i)->second;
        }
        gflow->addModifier(new BirthRate(gflow, birthRates));
      }
    }

    // --- Rules for destroying particles
    parser.getHeading_Optional("Destruction");
    hd = parser.first();
    if (hd) {
      // The subheads have no bodies, and are of the form: [type i] : [rate].
      // For now, we only allow reproduction to create the same type.
      std::map<int, RealType> death;
      for (auto s : hd->subHeads) {
        if (s->params.size()!=1) 
          throw BadStructure("Expected a rate (1 option). Found "+toStr(s->params.size())+".");
        death.insert(std::pair<int, RealType>(
          convert<int>(s->heading),
          convert<RealType>(s->params[0]->partA)
        ));
      }
      // Create a vector of birth rates - we must already know NTypes.
      if (!death.empty()) {
        vector<RealType> deathRates(NTypes, 0);
        for (int i=0; i<NTypes; ++i) {
          if (contains(death, i)) deathRates[i] = death.find(i)->second;
        }
        gflow->addModifier(new DeathRate(gflow, deathRates));
      }
    }

    parser.getHeading_Optional("Modifier");
    for (auto h : parser) add_modifier(h);

    // --- Get a maximum domain update delay time
    parser.getHeading_Optional("MaxDomainUpdateDelay");
    hd = parser.first();
    if (hd) {
      RealType max = -1;
      parser.extract_first_parameter(max);
      if (max>0) gflow->domain->setMaxUpdateDelay(max);
    }

    // --- Get a maximum allowed timestep
    parser.getHeading_Optional("MaxDT");
    hd = parser.first();
    if (hd) {
      RealType max = -1;
      parser.extract_first_parameter(max);
      if (max>0) gflow->integrator->setMaxDT(max);
    }

    // --- Get a maximum allowed timestep
    parser.getHeading_Optional("MinDT");
    hd = parser.first();
    if (hd) {
      RealType min = -1;
      parser.extract_first_parameter(min);
      if (min>0) gflow->integrator->setMinDT(min);
    }

    // Initialize domain
    gflow->domain->initialize();

    // --- Look for particle reconcilliation (should we remove overlapping particles?). Must do this after domain initialization.
    parser.getHeading_Optional("Reconcile");
    for (auto m : parser) {
      if (m->params.empty()) throw BadStructure("Need a remove option.");
        if (m->params[0]->partA=="Remove") {
          if (m->params[0]->partB.empty()) gflow->removeOverlapping(2.); // Remove particles overlapping by a large amount
          else gflow->removeOverlapping(convert<RealType>(m->params[0]->partB));
        }
        else throw BadStructure("Unrecognized remove option, [" + m->params[0]->partA + "].");
    }

    // Reconstruct domain
    gflow->domain->construct();
  }

  inline Integrator* FileParseCreator::choose_integrator(HeadNode *head) const {
    string token = head->params[0]->partA;
    if      (token=="VelocityVerlet")       return new VelocityVerlet(gflow);
    else if (token=="OverdampedIntegrator") return new OverdampedIntegrator(gflow);
    else if (token=="OverdampedLangevin")   return new OverdampedLangevinIntegrator(gflow);
    else if (token=="LangevinIntegrator") {
      LangevinIntegrator *integrator = new LangevinIntegrator(gflow);
      if (!head->subHeads.empty()) {
        // Gather options
        std::map<string, HeadNode*> options;
        for (auto h : head->subHeads) 
          options.insert(std::pair<string, HeadNode*>(h->heading, h));
        // Check for temperature
        auto it = options.find("Temperature");
        if (it!=options.end()) {
          integrator->setTemperature(convert<RealType>(it->second->params[0]->partA));
        }
        // Check for viscosity
        it = options.find("Viscosity");
        if (it!=options.end()) {
          integrator->setViscosity(convert<RealType>(it->second->params[0]->partA));
        }
      }
      return integrator;
    }
    else throw UnexpectedOption();
  }

  inline Interaction* FileParseCreator::choose_interaction(HeadNode *head) const {
    string token = head->params[2]->partA;
    if      (token=="HardSphere")   return new HardSphere(gflow);
    if      (token=="HardSphereGeneral") return new HardSphereGeneral(gflow);
    else if (token=="LennardJones") return new LennardJones(gflow);
    /*
    else if (token=="Coulomb")      return new CoulumbForce(gflow);
    else if (token=="Consumption")  return new Consumption(gflow);
    */
    else if (token=="None")         return nullptr;
    else throw UnexpectedOption();
  }

  inline BCFlag FileParseCreator::choose_bc(string& token) const {
    if      (token=="Wrap") return BCFlag::WRAP;
    else if (token=="Refl") return BCFlag::REFL;
    else if (token=="Repl") return BCFlag::REPL;
    else throw UnexpectedOption();
  }

  inline void FileParseCreator::add_modifier(HeadNode *head) const {
    string token = head->params[0]->partA;
    if (token=="WindTunnel")
      gflow->addModifier(new WindTunnel(gflow, convert<RealType>(head->params[0]->partB)));
    else throw UnexpectedOption();
  }

  inline void FileParseCreator::fillArea(HeadNode *head) const {
    // Check if head is good
    if (head==nullptr) return;
    // For local particle templates
    std::map<string, ParticleTemplate> particle_templates;

    // Create a parse helper
    ParseHelper parser(head);
    // Declare valid options
    parser.addValidSubheading("Bounds");
    parser.addValidSubheading("Template");
    parser.addValidSubheading("Number");
    parser.addValidSubheading("Velocity");
    parser.addValidSubheading("Relax");
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
    FillBounds *bnds;
    bnds = getFillBounds(parser.first());

    // --- Local Particle Template. Defines "types" of particles, e.g. radius distribution, density/mass, etc.
    parser.getHeading_Optional("Template");
    for (auto h : parser) getParticleTemplate(h, particle_templates);

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
      subParser.sortOptions();
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
    RealType *Vs = new RealType[sim_dimensions];
    int velocityOption = 0; // Normal velocities by default.
    hd = parser.first();
    if (hd) {
      string opt;
      parser.extract_first_parameter(opt, hd);
      if (opt=="Zero") velocityOption = 1; // Zero velocity
      else {
        parser.set_vector_argument(Vs, hd, sim_dimensions);
        velocityOption = 1;
      }
      // --- Other options go here
    }

    // --- Relaxation. How long to relax the particles
    parser.getHeading_Optional("Relax");
    RealType relax_length = -1;
    hd = parser.first();
    if (hd) parser.set_scalar_argument(relax_length, hd);

    // --- Check that we have defined a good area
    if (!usePhi && number<=0 && singleType)
      throw BadStructure("If using a single type, we need a nonzero number of particles.");

    // --- We have found all the options we need so far. Fill the area.
    GFlow filler(sim_dimensions);
    filler.setBounds(bnds->getBounds());
    // Default bounds are repulsive
    filler.setAllBCs(BCFlag::REPL);
    // If fill area has full bounds in any dimension, use the natural bounds
    for (int d=0; d<sim_dimensions; ++d) 
      if (
        filler.getBounds().min[d]==gflow->getBounds().min[d] 
        && filler.getBounds().max[d]==gflow->getBounds().max[d]
      ) filler.setBC(d, gflow->getBCs()[d]);

    // Delete filler's force master
    delete filler.forceMaster;
    // Set the new force master
    filler.forceMaster = gflow->forceMaster; // Make sure the particles treat each other in the same way

    // Get the simdata
    SimData *simData = filler.simData;

    // --- Attraction towards the center of the domain
    parser.getHeading_Optional("Attraction");
    hd = parser.first();
    if (hd) {
      RealType g;
      parser.set_scalar_argument(g, hd);
      filler.setAttraction(g);
    }
    
    // --- Fill with particles
    RealType *X = new RealType[sim_dimensions], *V = new RealType[sim_dimensions], sigma(0.), im(0.);
    int type(0);
    zeroVec(V, sim_dimensions);
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
      RealType vol = 0, Vol = bnds->vol();
      while (vol/Vol < phi) {
        // Select a position for the particle (random uniform)
        bnds->pick_position(X);
        // Choose a type of particle to create
        int pt = choice(global_generator);
        ParticleTemplate &particle_creator = particle_template_numbers.empty() ? particle_templates[0] : template_vector.at(pt);
        // Select other characteristics
        particle_creator.createParticle(X, sigma, im, type, i, sim_dimensions);
        // Add the particle
        simData->addParticle(X, V, sigma, im, type);
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
          bnds->pick_position(X);
          // Select other characteristics
          particle_creator.createParticle(X, sigma, im, type, i, sim_dimensions);
          // Add the particle
          simData->addParticle(X, V, sigma, im, type);
        }
      }
    }

    // Print status
    build_message += "From Fill Area: Done with initial particle assigmnemt.\n"; 

    // --- Relax the simulation
    if (relax_length<0) relax_length = 0.1; // Default relax length is 0.1
    if (relax_length>0) build_message += "From Fill Area: Relaxing for " + toStr(relax_length) + ".\n";
    hs_relax(&filler, relax_length); // To make sure particles don't stop on top of one another
    // relax(&filler, 0.15);

    // Select a velocity
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

    // --- Fill gflow with the particles
    copyVec(Vs, V, sim_dimensions);
    for (int i=0; i<simData->number(); ++i) {
      // Extract the particle properties
      copyVec(simData->X(i), X, sim_dimensions);
      int type = simData->Type(i);
      RealType sigma = simData->Sg(i);
      RealType im = simData->Im(i);
      if (type!=-1) {
        // Select the velocity for the final particle
        if (velocityOption==0) select_velocity(V, X, sigma, im, type);
        else if (velocityOption==1); // We already set V to be Vs
        else zeroVec(V, sim_dimensions);
        // Infinitely heavy objects do not move.
        if (im==0) zeroVec(V, sim_dimensions);
        gflow->simData->addParticle(X, V, sigma, im, type);
      }
    }
      
    // So we don't delete the force master when filler cleans up
    filler.forceMaster = nullptr; 

    // Clean up
    delete [] Vs;
    delete [] X;
    delete [] V;
    delete bnds;
  }

  inline FillBounds* FileParseCreator::getFillBounds(HeadNode *h) const {
    FillBounds *fbnds = nullptr;
    // Check if we should use the full simulation bounds - the bounds are then rectangular
    if (h->subHeads.empty() && h->params.size()>0 && h->params[0]->partA=="Full")
      fbnds = new RectangularBounds(simBounds, sim_dimensions);

    // Spherical bounds
    else if (h->params.size()>0 && h->params[0]->partA=="Sphere") {
      if (h->subHeads.size()!=2)
        throw BadStructure("Expected position and radius for spherical bounds.");
      if (h->subHeads[0]->params.size()!=sim_dimensions)
        throw BadStructure("Vector needs the correct number of dimensions.");
      SphericalBounds *sbnds = new SphericalBounds(sim_dimensions);
      // Get the center
      for (int d=0; d<sim_dimensions; ++d)
        sbnds->center[d] = convert<float>(h->subHeads[0]->params[d]->partA);
      // Get the radius
      sbnds->radius = convert<float>(h->subHeads[1]->params[0]->partA);
      fbnds = sbnds;
    }
    // Otherwise, rectangular bounds
    else {
      if (h->subHeads.size()!=sim_dimensions) 
        throw BadStructure("Expected "+toStr(sim_dimensions)+" arguments, found "+toStr(h->subHeads.size()));
      // Set bounds
      RectangularBounds *rbnds = new RectangularBounds(sim_dimensions);
      for (int d=0; d<h->subHeads.size(); ++d) {
        // Check for well formed options.
        if (h->subHeads[d]->params.size()!=2 || !h->subHeads[d]->subHeads.empty()) 
          throw BadStructure("Bounds need a min and a max, we found "+toStr(h->subHeads[d]->params.size())+" parameters.");
        // Extract the bounds.
        rbnds->min[d] = convert<float>( h->subHeads[d]->params[0]->partA );
        rbnds->max[d] = convert<float>( h->subHeads[d]->params[1]->partA );
      }
      fbnds = rbnds;
    }

    return fbnds;
  }

  inline void FileParseCreator::createParticle(HeadNode *head) const {
    // Check if head is good
    if (head==nullptr) return;

    // Create a parse helper
    ParseHelper parser(head);
    // Declare valid options
    parser.addValidSubheading("Position");
    parser.addValidSubheading("Velocity");
    parser.addValidSubheading("Sigma");
    parser.addValidSubheading("Type");
    parser.addValidSubheading("Mass");
    parser.addValidSubheading("Modifier");
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

    // --- Sort options
    std::multimap<string, HeadNode*> options;
    for (auto h : head->subHeads)
      options.insert(std::pair<string, HeadNode*>(h->heading, h));
    // Container to collect options in
    std::vector<HeadNode*> container;

    // --- Look for options
    RealType *X = new RealType[sim_dimensions], *V = new RealType[sim_dimensions], sigma, im;
    int type;

    // --- Look for position
    parser.getHeading_Necessary("Position", "Particle needs a position!");
    hd = parser.first();
    parser.set_vector_argument(X, hd, sim_dimensions);

    // --- Look for velocity
    parser.getHeading_Optional("Velocity");
    hd = parser.first();
    if (hd) parser.set_vector_argument(V, hd, sim_dimensions);
    else zeroVec(V, sim_dimensions);

    // --- Look for sigma
    parser.getHeading_Necessary("Sigma", "Particle needs a sigma!");
    hd = parser.first();
    parser.set_scalar_argument(sigma, hd);

    // --- Look for type
    parser.getHeading_Optional("Type");
    hd = parser.first();
    if (hd) parser.set_scalar_argument(type, hd);
    else type = 0;

    // --- Look for mass
    parser.getHeading_Optional("Mass");
    hd = parser.first();
    if (hd) {
      RealType d = 1.;
      if (hd->params[0]->partA=="inf") im = 0;
      else if (parser.extract_parameter(hd, "Density", d)) im = 1./(d*sphere_volume(sigma, sim_dimensions));
      else if (parser.extract_parameter(hd, "Mass", d)) im = 1./d;
      else throw BadStructure("Unrecognized option for mass in particle creation.");
    }
    else im = 1./sphere_volume(sigma, sim_dimensions); // Default density is 1

    // --- Look for modifiers
    parser.getHeading_Optional("Modifier");
    int g_id = gflow->simData->getNextGlobalID();
    for (auto m=parser.begin(); m!=parser.end(); ++m) {
      if (m.first_param()=="CV") gflow->addModifier(new ConstantVelocity(gflow, g_id, V));
      else if (m.first_param()=="CV-D") {
        RealType D = m.first_arg<RealType>();
        gflow->addModifier(new ConstantVelocityDistance(gflow, g_id, V, D));
      }
      else BadStructure("Unrecognized modifer option, ["+m.first_param()+"].");
    }

    // Add the particle to the system
    gflow->simData->addParticle(X, V, sigma, im, type);

    // Clean up
    delete [] X;
    delete [] V;
  }

  inline void FileParseCreator::getAllMatches(string heading, vector<HeadNode*>& container, std::multimap<string, HeadNode*>& options) const {
    // Clear container in case there is stuff left over in it.
    container.clear();
    // Look for options
    for (auto it=options.find(heading); it!=options.end(); ++it) {
      // If the head has the propper heading, store it.
      if (it->first==heading) container.push_back(it->second);
    }
  }

  //! Valid headers for particle templates are:
  //! Radius: [Random | Uniform | Normal | #]
  //! Mass:   [Density=# | #]
  //! Type:   [Random | #]
  inline void FileParseCreator::getParticleTemplate(HeadNode *head, std::map<string, ParticleTemplate>& particle_templates) const {
    // Create a particle template to set
    ParticleTemplate p_template;

    // --- Sort options
    std::multimap<string, HeadNode*> options;
    for (auto h : head->subHeads)
      options.insert(std::pair<string, HeadNode*>(h->heading, h));
    // For collecting options
    std::vector<HeadNode*> container;
    string type;

    // --- Look for options
    
    // --- Look for Radius option
    getAllMatches("Radius", container, options);
    if (container.empty()) 
      throw BadStructure("We need some radius information!");
    if (container.size()>1)
      build_message += "We only need one radius information block. Ignoring all but the first instance.\n";
    // Get the random engine
    p_template.radius_engine = getRandomEngine(container[0], type);
    p_template.radius_string = type;
    
    // --- Look for Mass option
    getAllMatches("Mass", container, options);
    if (container.empty())
      throw BadStructure("We need some mass information!");
    if (container.size()>1)
      build_message += "We only need one mass information block. Ignoring all but the first instance.\n";
    // Get Density
    if (container[0]->params[0]->partA=="Density") {
      p_template.mass_engine = new DeterministicEngine(convert<RealType>(container[0]->params[0]->partB));
      p_template.mass_string = "Density";
    }
    // Otherwise, get a number
    else {
      p_template.mass_engine = new DeterministicEngine(convert<RealType>(container[0]->params[0]->partB));
      p_template.mass_string = "Mass";
    }
    //p_template.mass_engine = getRandomEngine(container[0], type);
    //p_template.mass_string = type;

    // --- Look for Type option
    getAllMatches("Type", container, options);
    if (container.empty())
      throw BadStructure("We need some type information!");
    if (container.size()>1)
      build_message += "We only need one type information block. Ignoring all but the first instance.\n";
    p_template.type_engine = getRandomEngine(container[0], type);
    p_template.type_string = type;

    // Add to particle templates
    particle_templates.insert(std::pair<string, ParticleTemplate>(head->params[0]->partA, p_template));
  }

  inline void FileParseCreator::makeRandomForces() {
    // Assign random interactions, either LennardJones or HardSphere (for now), and with equal probability (for now)
    Interaction *hardSphere   = new HardSphere(gflow);
    Interaction *lennardJones = new LennardJones(gflow);
    // Assign random (but symmetric) interactions
    for (int i=0; i<NTypes; ++i) {
      // Self interaction
      if (drand48()>0.5) gflow->forceMaster->setInteraction(i, i, hardSphere);
      else gflow->forceMaster->setInteraction(i, i, lennardJones);

      for (int j=i+1; j<NTypes; ++j) {
        if (drand48()>0.5) {
          gflow->forceMaster->setInteraction(i, j, hardSphere);
          gflow->forceMaster->setInteraction(j, i, hardSphere);
        }
        else {
          gflow->forceMaster->setInteraction(i, j, lennardJones);
          gflow->forceMaster->setInteraction(j, i, lennardJones);
        }
      }
    }
  }

  inline RandomEngine* FileParseCreator::getRandomEngine(HeadNode *h, string &type) const {
    // Clear type string
    type.clear();
    // Check structure
    if (h->params.size()!=1)
      throw BadStructure("Random engine needs one parameter, found "+toStr(h->params.size())+".");
    else if (!h->params[0]->partB.empty()) {
      // A type is specified
      type = h->params[0]->partA;
      // We expect a number in partB
      return new DeterministicEngine(convert<double>(h->params[0]->partB));
    }
    // If there is no body
    else if (h->subHeads.empty()) {
      string token = h->params[0]->partA;
      // We expect either a number in partA, or a string (e.g. "Inf").
      if (isdigit(token.at(0))) {
        return new DeterministicEngine(convert<double>(token));
      }
      else if (h->heading=="Type" && token=="Equiprobable") {
        type = token;
        // Only for type - Create each type with equal probability
        RealType p = 1./NTypes;
        vector<double> probabilities, values;
        for (int i=0; i<NTypes; ++i) {
          probabilities.push_back(p);
          values.push_back(i);
        }
        return new DiscreteRandomEngine(probabilities, values);
      }
      else {
        type = token;
        // The engine will not matter
        return new DeterministicEngine(0);
      }
    }

    // Parameter string
    string param = h->params[0]->partA;
    // Check for options:
    if (param=="Random") {
      // Discrete Random, with specified values and probabilities
      vector<double> probabilities, values;
      for (auto sh : h->subHeads) {
        if (sh->heading!="P")
          throw BadStructure("Discrete probability should be indicated with a 'P.'");
        if (sh->params.size()!=2)
          throw BadStructure("Expect value and probability, encountered "+toStr(sh->params.size())+" options.");
        if (!sh->params[0]->partB.empty() || !sh->params[1]->partB.empty())
          throw BadStructure("Expected one part parameters in discrete probability specification.");
        // Store the value and its probability
        values.push_back(convert<double>(sh->params[0]->partA)); 
        probabilities.push_back(convert<double>(sh->params[1]->partA));
      }
      return new DiscreteRandomEngine(probabilities, values);
    }
    else if (param=="Uniform") {
      if (h->subHeads.size()!=2) 
        throw BadStructure("Uniform distribution needs a min and a max. Found "+toStr(h->subHeads.size())+" parameters.");
      HeadNode *lwr = h->subHeads[0], *upr = h->subHeads[1];
      // Swap if necessary
      if (lwr->heading!="Min") std::swap(lwr, upr);
      // Check structure
      if (lwr->heading!="Min" || upr->heading!="Max") 
        throw BadStructure("Headings incorrect for uniform distribution. Need [Min] and [Max]");
      if (lwr->params.size()!=1 || upr->params.size()!=1) 
        throw BadStructure("More than one parameter for uniform distribution.");
      // Set the uniform random engine
      return new UniformRandomEngine(
        convert<double>(lwr->params[0]->partA), convert<double>(upr->params[0]->partA)
      );
    }
    else if (param=="Normal") {
      if (h->subHeads.size()!=2) 
        throw BadStructure("Normal distribution needs average and variance. Found "+toStr(h->subHeads.size())+" parameters.");
      HeadNode *ave = h->subHeads[0], *var = h->subHeads[1];
      // Swap if necessary
      if (ave->heading!="Ave") std::swap(ave, var);
      // Check structure
      if (ave->heading!="Ave" || var->heading!="Var")
        throw BadStructure("Headings incorrect for normal distribution.");
      if (ave->params.size()!=1 || var->params.size()!=1) 
        throw BadStructure("More than one parameter for normal distribution.");
      // Set the uniform random engine
      return new NormalRandomEngine(
        convert<double>(ave->params[0]->partA), convert<double>(var->params[0]->partA)
      );
    }
    else throw BadStructure("Unrecognized choice for a random engine.");
    // We should never reach here
    throw BadStructure("An unreachable part of code was reached!");
    // Token return
    return nullptr;
  }

  inline string FileParseCreator::copyFile() const {
    string s;
    char c;
    ifstream fin(configFile);
    if (fin.fail()) return "";
    fin.get(c);
    while (!fin.eof()) {
      s.push_back(c);
      fin.get(c);
    }
    return s;
  }

}
