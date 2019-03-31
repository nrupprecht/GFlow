#include "fileparsecreator.hpp"
// Other files
#include "../utility/printingutility.hpp"
#include "../allmodifiers.hpp"
#include "../interactions/interaction-choice.hpp"
#include "../allareacreators.hpp"

namespace GFlowSimulation {

  FileParseCreator::FileParseCreator(int argc, char **argv) : Creator(argc, argv), configFile(""), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    // Seed generators here
    seedGenerator(seed);
    generator = std::mt19937(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  FileParseCreator::FileParseCreator(ArgParse *p) : Creator(p), configFile(""), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    // Seed generators here
    seedGenerator(seed);
    generator = std::mt19937(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  FileParseCreator::FileParseCreator(ArgParse *p, string f) : Creator(p), configFile(f), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    // Seed generators here
    seedGenerator(seed);
    generator = std::mt19937(seed);
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
    catch (UnexpectedOption uo) {
      cout << "Caught Unexpected Option error: " << uo.message << endl;
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
    parser.addValidSubheading("Var");
    parser.addValidSubheading("Dimensions");
    parser.addValidSubheading("Seed");
    parser.addValidSubheading("Bounds");
    parser.addValidSubheading("Boundary");
    parser.addValidSubheading("Gravity");
    parser.addValidSubheading("Attraction");
    parser.addValidSubheading("Integrator");
    parser.addValidSubheading(Types_Token);
    parser.addValidSubheading(Interactions_Token);
    parser.addValidSubheading("Template");
    parser.addValidSubheading("Fill");
    parser.addValidSubheading("Particle");
    parser.addValidSubheading("Creation");
    parser.addValidSubheading("Destruction");
    parser.addValidSubheading("Modifier");
    parser.addValidSubheading("MaxDomainUpdateDelay");
    parser.addValidSubheading("MaxDT");
    parser.addValidSubheading("MinDT");
    parser.addValidSubheading("Relax");
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

    // --- Look for variables
    parser.getHeading_Optional("Var");
    for (auto v=parser.begin(); v!=parser.end(); ++v) {
      string name, value;
      name  = v.first_param();
      value = v.first_arg<string>();
      // Check for command line argument
      parserPtr->get(name, value);
      // Add to variables
      variables.insert(pair<string, string>(name, value));
    }
    // Set the variables
    parser.set_variables(variables);

    // --- Look for dimensions
    parser.getHeading_Optional("Dimensions");
    hd = parser.first();
    sim_dimensions = 2; // Default value
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
    subParser.set_variables(variables);
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
      if (parser.extract_first_parameter(opt)) gflow->setBC(d, choose_bc(opt));
      else {
        ParseHelper subParser(hd);
        subParser.set_variables(variables);
        // Record template name, number
        int d=0;
        for (auto m=subParser.begin() ; m!=subParser.end(); ++m, ++d) {
          opt = m.convert_param<string>();
          gflow->setBC(d, choose_bc(opt));
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
    parser.getHeading_Optional("Fill");
    FillAreaCreator fillarea;
    PolymerCreator pc;
    for (auto h : parser) {
      if (h->params[0]->partA=="Area") fillarea.createArea(h, gflow, variables, particle_fixers);
      else if (h->params[0]->partA=="Polymer") pc.createArea(hd, gflow, variables, particle_fixers);
    }

    // --- Create a single particle
    parser.getHeading_Optional("Particle");
    for (auto h : parser) createParticle(h);

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

    // --- Relax the simulation
    parser.getHeading_Optional("Relax");
    hd = parser.first();
    // Default relax time is 0.5
    RealType relax_time = 0.5;
    if (hd) {
      parser.extract_first_parameter(relax_time);
    }
    if (relax_time>0) hs_relax(gflow, relax_time);

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

    // Fix velocities
    fix_particle_velocities(gflow->simData);
  }

  inline Integrator* FileParseCreator::choose_integrator(HeadNode *head) const {
    string token = head->params[0]->partA;
    Integrator *integrator = nullptr;
    if      (token=="VelocityVerlet")       integrator = new VelocityVerlet(gflow);
    else if (token=="OverdampedIntegrator") integrator = new OverdampedIntegrator(gflow);
    else if (token=="OverdampedLangevin")   integrator = new OverdampedLangevinIntegrator(gflow);
    else if (token=="LangevinIntegrator")    integrator = new LangevinIntegrator(gflow);
    else throw UnexpectedOption("Integrator choice was ["+token+"].");
    // Temperature and viscosity for LangevinType integrators
    LangevinTypeIntegrator* lti = dynamic_cast<LangevinTypeIntegrator*>(integrator);
    if (lti) {
      if (!head->subHeads.empty()) {
        // Gather options
        std::map<string, HeadNode*> options;
        for (auto h : head->subHeads) 
          options.insert(std::pair<string, HeadNode*>(h->heading, h));
        // Check for temperature
        auto it = options.find("Temperature");
        if (it!=options.end()) {
          lti->setTemperature(convert<RealType>(it->second->params[0]->partA));
        }
        // Check for viscosity
        it = options.find("Viscosity");
        if (it!=options.end()) {
          lti->setViscosity(convert<RealType>(it->second->params[0]->partA));
        }
      }
    }
    // Look for general options
    ParseHelper parser(head);
    parser.set_variables(variables);
    parser.addValidSubheading("Delay");
    parser.addValidSubheading("MaxDT");
    parser.addValidSubheading("MinDT");
    parser.addValidSubheading("Adjust");
    parser.addValidSubheading("UseV");
    parser.addValidSubheading("UseA");
    parser.sortOptions();
    HeadNode *hd = nullptr;

    parser.getHeading_Optional("Delay");
    hd = parser.first();
    if (hd) {
      int steps = 0;
      parser.extract_first_parameter(steps, hd);
      integrator->setStepDelay(steps);
    }

    parser.getHeading_Optional("MaxDT");
    hd = parser.first();
    if (hd) {
      RealType max_dt;
      parser.extract_first_parameter(max_dt, hd);
      if (max_dt>0) integrator->setMaxDT(max_dt);
    }

    parser.getHeading_Optional("MinDT");
    hd = parser.first();
    if (hd) {
      RealType min_dt;
      parser.extract_first_parameter(min_dt, hd); 
      if (min_dt>0) integrator->setMinDT(min_dt);
    }

    parser.getHeading_Optional("Adjust");
    hd = parser.first();
    if (hd) {
      int adj;
      parser.extract_first_parameter(adj, hd); 
      if (adj>=0) integrator->setAdjustDT(adj);
    }

    parser.getHeading_Optional("UseV");
    hd = parser.first();
    if (hd) {
      int u;
      parser.extract_first_parameter(u, hd); 
      if (u>=0) integrator->setUseV(u);
    }

    parser.getHeading_Optional("UseA");
    hd = parser.first();
    if (hd) {
      int u;
      parser.extract_first_parameter(u, hd); 
      if (u>=0) integrator->setUseA(u);
    }

    // Return
    return integrator;
  }

  inline Interaction* FileParseCreator::choose_interaction(HeadNode *head) const {
    string token = head->params[2]->partA;
    return InteractionChoice::choose(gflow, token, sim_dimensions);
  }

  inline BCFlag FileParseCreator::choose_bc(string& token) const {
    if      (token=="Wrap")     return BCFlag::WRAP;
    else if (token=="Reflect")  return BCFlag::REFL;
    else if (token=="Repulse") return BCFlag::REPL;
    else throw UnexpectedOption("Boundary condition choice was ["+token+"].");
  }

  inline void FileParseCreator::add_modifier(HeadNode *head) const {
    ParseHelper parser(nullptr);
    parser.set_variables(variables);
    string token = head->params[0]->partA;
    if (token=="WindTunnel")
      gflow->addModifier(new WindTunnel(gflow, parser.value<RealType>(head->params[0]->partB)));
    else throw UnexpectedOption("Modifier choice was ["+token+"].");
  }

  inline void FileParseCreator::createParticle(HeadNode *head) const {
    // Check if head is good
    if (head==nullptr) return;

    // Create a parse helper
    ParseHelper parser(head);
    parser.set_variables(variables);
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

    ParseHelper parser(head);
    parser.set_variables(variables);
    // Declare valid options
    parser.addValidSubheading("Radius");
    parser.addValidSubheading("Mass");
    parser.addValidSubheading("Type");
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
    
    // --- Look for options
    string option, type;
    
    // --- Look for Radius option
    parser.getHeading_Necessary("Radius");
    hd = parser.first();
    p_template.radius_engine = getRandomEngine(hd, type);
    p_template.radius_string = type;

    // --- Look for Mass option
    parser.getHeading_Necessary("Mass");
    hd = parser.first();
    parser.extract_first_parameter(option, hd);
    RealType m;
    if (option=="Density" && parser.extract_first_arg(m)) {
      p_template.mass_engine = new DeterministicEngine(m);
      p_template.mass_string = "Density";
    }
    else if (parser.extract_first_arg(m)) {
      p_template.mass_engine = new DeterministicEngine(m);
      p_template.mass_string = "Mass";
    }

    // --- Look for Type option
    parser.getHeading_Necessary("Type");
    hd = parser.first();
    p_template.type_engine = getRandomEngine(hd, type);
    p_template.type_string = type;

    // Add to particle templates
    particle_templates.insert(std::pair<string, ParticleTemplate>(head->params[0]->partA, p_template));
  }

  inline void FileParseCreator::makeRandomForces() {
    // Assign random interactions, either LennardJones or HardSphere (for now), and with equal probability (for now)
    Interaction *hs   = InteractionChoice::choose(gflow, InteractionChoice::HardSphereToken, sim_dimensions); 
    Interaction *lj = InteractionChoice::choose(gflow, InteractionChoice::LennardJonesToken, sim_dimensions);
    // Assign random (but symmetric) interactions
    for (int i=0; i<NTypes; ++i) {
      // Self interaction
      if (drand48()>0.5) gflow->forceMaster->setInteraction(i, i, hs);
      else gflow->forceMaster->setInteraction(i, i, lj);

      for (int j=i+1; j<NTypes; ++j) {
        if (drand48()>0.5) {
          gflow->forceMaster->setInteraction(i, j, hs);
          gflow->forceMaster->setInteraction(j, i, hs);
        }
        else {
          gflow->forceMaster->setInteraction(i, j, lj);
          gflow->forceMaster->setInteraction(j, i, lj);
        }
      }
    }
  }

  inline RandomEngine* FileParseCreator::getRandomEngine(HeadNode *h, string &type) const {
    // Clear type string
    type.clear();

    ParseHelper parser(h);
    parser.set_variables(variables);

    // Check structure
    if (h->params.size()!=1)
      throw BadStructure("Random engine needs one parameter, found "+toStr(h->params.size())+".");
    else if (!h->params[0]->partB.empty()) {
      // A type is specified
      type = h->params[0]->partA;
      // We expect a number in partB
      return new DeterministicEngine(parser.value<double>(h->params[0]->partB));
    }
    // If there is no body
    else if (h->subHeads.empty()) {
      string token = h->params[0]->partA;
      // We expect either a number in partA, or a string (e.g. "Inf").
      if (isdigit(token.at(0))) {
        return new DeterministicEngine(parser.value<double>(token));
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
        values.push_back(parser.value<double>(sh->params[0]->partA)); 
        probabilities.push_back(parser.value<double>(sh->params[1]->partA));
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
        parser.value<double>(lwr->params[0]->partA), parser.value<double>(upr->params[0]->partA)
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
        parser.value<double>(ave->params[0]->partA), parser.value<double>(var->params[0]->partA)
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
