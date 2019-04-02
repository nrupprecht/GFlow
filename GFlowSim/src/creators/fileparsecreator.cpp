#include "fileparsecreator.hpp"
// Other files
#include "../utility/printingutility.hpp"
#include "../allmodifiers.hpp"
#include "../interactions/interaction-choice.hpp"
#include "../allareacreators.hpp"

#include "../utility/treeparser.hpp"

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
    bool parse_trace = false;
    // Get command line arguments
    if (parserPtr) {
      parserPtr->get("parse-trace", parse_trace);
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
    // Create from the options
    build_message += "Starting simulation setup... ";
    try {
      createFromOptions(root);
    }
    catch (BadStructure bs) {
      cout << "Caught Bad Structure error: Message: [" << bs.message << "]\n";
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

    // If requested, print out the parse trace.
    if (parse_trace) {
      // Print the parse tree to a file so we can debug
      ofstream fout("ParseTrace.txt");
      if (fout.fail());
      else {
        fout << parse_message;
        fout.close();
      }
    }

    // Clean up
    delete root;
    
    // Timing
    auto end_time = current_time();
    gflow->dataMaster->setInitializationTime(time_span(end_time, start_time));
    // Give the setup file to the datamaster.
    gflow->dataMaster->giveFile("setup.txt", file);

    // Clean up and return
    GFlow *ret = gflow;
    gflow = nullptr;
    return ret;
  }

  void FileParseCreator::setVariable(const string& name, const string& value) {
    // Insert a variable if it does not already exist.
    if (variables.find(name)==variables.end())
      variables.insert(pair<string, string>(name, value));
  }

  inline void FileParseCreator::createFromOptions(HeadNode *head) {
    // Create a parser
    TreeParser parser(head);
    // Declare valid options.
    parser.addHeadingNecessary("Bounds", "We need bounds!");
    parser.addHeadingNecessary(Interactions_Token, "We need interaction information!");
    parser.addHeadingOptional("Var");
    parser.addHeadingOptional("Dimensions");
    parser.addHeadingOptional("Seed");
    parser.addHeadingOptional("Boundary");
    parser.addHeadingOptional("Gravity");
    parser.addHeadingOptional("Attraction");
    parser.addHeadingOptional("Integrator");
    parser.addHeadingOptional(Types_Token);
    parser.addHeadingOptional("Template");
    parser.addHeadingOptional("Fill");
    parser.addHeadingOptional("Particle");
    parser.addHeadingOptional("Modifier");
    parser.addHeadingOptional("Relax");
    parser.addHeadingOptional("Reconcile");
    // Check headings for validity.
    parser.check();

    // --- Look for variables
    if (parser.begin("Var")) {
      string name, value;
      do {
        // Get the name of the variable, and its value.
        parser.argvalName(name, value);
        // Check for command line argument
        parserPtr->get(name, value);
        // Add to variables if a variable does not already exist.
        setVariable(name, value);
      } while (parser.next());
      // Set the variables
      parser.set_variables(variables);
    }

    // --- Look for dimensions
    sim_dimensions = 2;
    parser.firstArg("Dimensions", sim_dimensions);
    if (sim_dimensions<=0) throw BadDimension();

    // Set up GFlow
    if (gflow) delete gflow;
    gflow = new GFlow(sim_dimensions);

    // --- Look for seed information for random generators
    if (parser.firstArg("Seed", seed));
    else seed = std::chrono::system_clock::now().time_since_epoch().count();
    // Set up random number generators
    seedGenerator(seed);
    srand48(seed);
    normal_dist = std::normal_distribution<RealType>(0., 1.);

    // --- Look for bounds
    parser.focus("Bounds"); // This will be true, since Bounds is a required heading.
    // Parse to find the bounds
    Region *region = ParseConstructor::parse_region(parser.getNode(), variables, gflow);
    gflow->setBounds(region->get_bounding_box());
    delete region;
    // Return to original level
    parser.up();

    // --- Look for boundary conditions
    if (parser.focus("Boundary")) {
      // Set all the boundary conditions to be the same.
      if (parser.argName()!="") gflow->setAllBCs(choose_bc(parser.argName()));
      // Each boundary condition is different.
      else if (parser.body_size()==sim_dimensions) {
        // Get all the boundary conditions
        int d = 0;
        for (parser.begin(); parser.next(); ++d) {
          gflow->setBC(d, choose_bc(parser.argName()));
        }
      }
      else throw BadDimension();
      // Return to original level
      parser.up();
    }
    else gflow->setAllBCs(BCFlag::WRAP);

    // --- Constant, uniform acceleration.
    Vec g(0);
    parser.argVec("Gravity");
    if (g.size()!=0) gflow->addModifier(new ConstantAcceleration(gflow, g.data));

    // --- Attraction towards the center of the domain
    RealType att = 0;
    if (parser.firstArg("Attraction", att)) gflow->setAttraction(att);

    // --- Look for integrator
    if (parser.focus("Integrator")) {
      gflow->integrator = choose_integrator(parser.getNode());
      // Return to original level
      parser.up();
    }
    else gflow->integrator = new VelocityVerlet(gflow);

    // --- Look for number of particle types
    if (parser.firstArg(Types_Token, NTypes));
    else NTypes = 1;
    gflow->forceMaster->setNTypes(NTypes);

    // --- Look for interactions
    parser.focus(Interactions_Token);
    // An interaction grid is not specified.
    if (parser.body_size()==0) {
      // If all types interact with the same interaction, or some interaction command is given.
      if (parser.argName()=="Random") { 
        makeRandomForces();
      }
      // Expect the force type. Every particle type interacts via the same interaction.
      else {
        gflow->forceMaster->setInteraction(InteractionChoice::choose(gflow, parser.argName(), sim_dimensions));
      }
    }
    // The interaction grid is specified.
    else { 
      std::map<string, Interaction*> interactions;
      string token;
      int t1, t2;
      parser.begin();
      // Go through each entry in the force table (each head is an entry).
      do {
        // CASE: [ : type1, type2, interaction-type ]
        if (parser.args_size()==3) {
          // Get the interaction token
          token = parser.argName(2);
          // Get the particle types
          t1 = parser.arg_cast<int>(0);
          t2 = parser.arg_cast<int>(1);
          // If the interaction has not occured yet, create one.
          if (interactions.find(token)==interactions.end()) 
            interactions.insert(std::pair<string, Interaction*>(token, InteractionChoice::choose(gflow, token, sim_dimensions)));
          // Set the interparticle interaction
          gflow->forceMaster->setInteraction(t1, t2, interactions.find(token)->second);
        }
        // NOT A VALID CASE
        else throw BadStructure("Force grid: unrecognized structure.");
      } while (parser.next());
      // We have scanned all the entries.
    }
    // Return to the original level.
    parser.up();

    // --- Global Particle Template. Defines "types" of particles, e.g. radius distribution, density/mass, etc.
    if (parser.begin("Template")) {
      do {
        getParticleTemplate(parser.getNode(), global_templates);
      } while (parser.next());
    }

    // --- Fill an area with particles
    if (parser.begin("Fill")) {
      FillAreaCreator fa;
      PolymerCreator pc;
      do {
        if      (parser.argName()=="Area")    fa.createArea(parser.getNode(), gflow, variables, particle_fixers);
        else if (parser.argName()=="Polymer") pc.createArea(parser.getNode(), gflow, variables, particle_fixers);
        else throw BadStructure("Unrecognized fill option: [" + parser.argName() + "].");
      } while (parser.next());
    }

    // --- Create single particles
    if (parser.begin("Particle")) {
      do {
        createParticle(parser.getNode());
      } while (parser.next());
    }

    // --- Create single particles
    if (parser.begin("Modifier")) {
      do {
        add_modifier(parser.getNode());
      } while (parser.next());
    }

    // Initialize domain
    gflow->domain->initialize();

    // --- Relax the simulation
    RealType relax_time = 0.5;
    parser.firstArg("Relax", relax_time);
    if (relax_time>0) hs_relax(gflow, relax_time);

    // --- Reconcile overlaps
    if (parser.begin("Reconcile")) {
      do {
        if (parser.argName()=="Remove") {
          // Default remove overlap cutoff.
          RealType rem = 2.;
          // Remove overlap cutoff
          parser.val(rem);
          // Do the removal
          gflow->removeOverlapping(rem);
        }
        else throw BadStructure("Unrecognized remove option, [" + parser.argName() + "].");
      } while (parser.next());
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

  inline BCFlag FileParseCreator::choose_bc(const string& token) const {
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
