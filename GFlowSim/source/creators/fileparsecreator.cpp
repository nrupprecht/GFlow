#include <creators/fileparsecreator.hpp>
// Other files
#include <utility/printingutility.hpp>
#include <allmodifiers.hpp>
#include <interactions/interaction-choice.hpp>
#include <allareacreators.hpp>
#include <utility>

using namespace GFlowSimulation;

FileParseCreator::FileParseCreator(int argc, char **argv)
    : Creator(argc, argv), configFile(""), gflow(nullptr) {
  seed = std::chrono::system_clock::now().time_since_epoch().count();
  // Seed generators here
  seedGenerator(seed);
  generator = std::mt19937(seed);
  normal_dist = std::normal_distribution<RealType>(0., 1.);
}

FileParseCreator::FileParseCreator(ArgParse *p)
    : Creator(p), configFile(""), gflow(nullptr) {
  seed = std::chrono::system_clock::now().time_since_epoch().count();
  // Seed generators here
  seedGenerator(seed);
  generator = std::mt19937(seed);
  normal_dist = std::normal_distribution<RealType>(0., 1.);
}

FileParseCreator::FileParseCreator(ArgParse *p, string f)
    : Creator(p), configFile(std::move(f)), gflow(nullptr) {
  seed = std::chrono::system_clock::now().time_since_epoch().count();
  // Seed generators here
  seedGenerator(seed);
  generator = std::mt19937(seed);
  normal_dist = std::normal_distribution<RealType>(0., 1.);
  // Create the parse tree from the config file
  parseFile();
}

FileParseCreator::~FileParseCreator() {
  delete root;
  root = nullptr;
}

GFlow *FileParseCreator::createSimulation() {
  // Timing
  auto start_time = current_time();

  // Command line arguments
  bool parse_trace = false;
  // Get command line arguments
  if (parserPtr) {
    parserPtr->get("parse-trace", parse_trace);
  }

  if (root == nullptr) {
    if (!parseFile()) {
      return nullptr;
    }
  }

  // If requested, print out the parse trace.
  if (parse_trace) {
    // Print the parse tree to a file so we can debug
    ofstream fout("ParseTrace.txt");
    if (fout.fail()) {
    }
    else {
      fout << parse_message;
      fout.close();
    }
  }

  // Create the scenario from the options
  delete gflow;
  // Create from the options
  build_message += "Starting simulation setup... ";
  try {
    createFromOptions(root);
  }
  catch (BadStructure &bs) {
    cout << "--> Caught Bad Structure error: Message: [" << bs.message << "]\n";
    cout << "Build Message:\n" << build_message << endl;
    // Print the parse tree to a file so we can debug
    ofstream fout("ParseTrace.txt");
    if (fout.fail()) {
    }
    else {
      fout << parse_message;
      fout.close();
    }
    throw;
  }
  catch (UnexpectedOption &uo) {
    cout << "--> Caught Unexpected Option error: " << uo.message << endl;
    cout << "Build Message:\n" << build_message << endl;
    // Print the parse tree to a file so we can debug
    ofstream fout("ParseTrace.txt");
    if (fout.fail()) {
    }
    else {
      fout << parse_message;
      fout.close();
    }
    throw;
  }
  build_message += "Done.\n";

  //! \brief Collect the global ids.
  Creator::correct_global_ids(gflow);

  // Timing
  auto end_time = current_time();
  gflow->dataMaster->setInitializationTime(time_span(end_time, start_time));

  // Give the setup file to the datamaster.
  gflow->dataMaster->giveFile("setup.txt", setup_file);

  // Clean up and return
  GFlow *ret = gflow;
  gflow = nullptr;
  return ret;
}

void FileParseCreator::setVariable(const string &name, const string &value, bool overwrite) {
  // Try to find the variable
  auto var = variables.find(name);
  // Insert a variable if it does not already exist.
  if (var == variables.end()) {
    variables.insert(pair<string, string>(name, value));
    // Else, if overwrite is true, overwrite the existing variable
  }
  else if (overwrite) {
    var->second = value;
  }
}

inline bool FileParseCreator::parseFile() {
  // There needs to be a config file.
  if (configFile.empty()) {
    return false;
  }
  // We treat the file as one giant body. Get that body.
  build_message = "Starting file parse... ";
  // Create a file parser
  FileParse parser;
  // Clean up old tree
  delete root;
  // Try to parse the file.
  try {
    root = parser.parseFile(configFile); // A file parse class does the parsing
  }
  catch (FileParse::UnexpectedToken &ut) {
    cout << "Caught unexpected token error from file parsing. Message: " << ut.message << endl;
    cout << " Trace:\n";
    cout << parser.getMessage();
    throw;
  }
  build_message += "Done.\n";
  // Get the message from the parser.
  parse_message = parser.getMessage();
  // Save the file we just loaded
  setup_file = copyFile();
  // Return success.
  return true;
}

inline void FileParseCreator::createFromOptions(HeadNode *head) {
  // Create a parser
  TreeParser parser(head);
  // Declare valid options.
  parser.addHeadingNecessary("Bounds", "We need bounds!");
  parser.addHeadingNecessary(Interactions_Token, "We need interaction information!");
  parser.addHeadingOptional("Var");
  parser.addHeadingOptional("Value");
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
  parser.addHeadingOptional("HSRelax");
  parser.addHeadingOptional("Relax");
  parser.addHeadingOptional("Reconcile");
  parser.addHeadingOptional("SkinDepth");
  // Check headings for validity.
  parser.check();

  // Clear particle fixers
  clear_particle_fixers();

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

  // --- Look for values
  if (parser.begin("Value")) {
    string name, value;
    do {
      // Get the name of the variable, and its value.
      parser.argvalName(name, value);
      // Check for command line argument
      parserPtr->get(name, value);
      // If there is no value, something is wrong
      if (value.empty()) {
        throw BadStructure("Value definition should not be empty for value: [" + name + "].");
      }
      // Evaluate the value
      RealType val = Eval::evaluate(value, variables);
      // Add to variables if a variable does not already exist.
      setVariable(name, toStr(val));
      // Set the variables each time, so later values can depend on earlier values
      parser.set_variables(variables);
    } while (parser.next());

  }

  // --- Look for dimensions
  sim_dimensions = 2;
  parser.firstArg("Dimensions", sim_dimensions);
  if (sim_dimensions <= 0) {
    throw BadDimension();
  }

  // Set up GFlow
  delete gflow;
  gflow = new GFlow(sim_dimensions);

  // --- Look for seed information for random generators
  if (parser.firstArg("Seed", seed)) {
  }
  else {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
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
    if (!parser.argName().empty()) {
      gflow->setAllBCs(choose_bc(parser.argName()));
      // Each boundary condition is different.
    }
    else if (parser.body_size() == sim_dimensions) {
      // Get all the boundary conditions
      int d = 0;
      parser.begin();
      do {
        gflow->setBC(d, choose_bc(parser.argName()));
        ++d;
      } while (parser.next());
    }
    else {
      throw BadDimension();
    }
    // Return to original level
    parser.up();
  }
  else {
    gflow->setAllBCs(BCFlag::WRAP);
  }

  // --- Constant, uniform acceleration.
  Vec g = parser.argVec("Gravity");
  if (g.size() != 0) {
    gflow->addModifier(make_shared<ConstantAcceleration>(gflow, g.data));
  }

  // --- Attraction towards the center of the bounds
  RealType att = 0;
  if (parser.firstArg("Attraction", att)) {
    gflow->setAttraction(att);
  }

  // --- Look for integrator
  if (parser.focus("Integrator")) {
    gflow->integrator = choose_integrator(parser.getNode());
    // Return to original level
    parser.up();
  }
  else {
    gflow->integrator = choose_velocity_verlet(gflow, sim_dimensions);
  }

  // --- Look for number of particle types
  if (parser.firstArg(Types_Token, NTypes)) {
  }
  else {
    NTypes = 1;
  }
  gflow->forceMaster->setNTypes(NTypes);

  // --- Look for interactions
  parser.focus(Interactions_Token);
  // An interaction grid is not specified.
  if (parser.body_size() == 0) {
    // If all types interact with the same interaction, or some interaction command is given.
    if (parser.argName() == "Random") {
      makeRandomForces();
      // Expect the force type. Every particle type interacts via the same interaction.
    }
    else {
      // This allows the force string to be a variable.
      string token = parser.get_variable_string(parser.argName());
      if (token.empty()) {
        token = parser.argName();
      }
      // Set all forces to be the chosen force.
      gflow->forceMaster->setInteraction(ParseConstructor::getInteraction(parser.getNode(), variables, token, gflow));
    }
  }
    // Expect the force type, with a body specifying parameters. Every particle type interacts via the same interaction.
  else if (parser.args_size() == 1) {
    // This allows the force string to be a variable.
    string token = parser.get_variable_string(parser.argName());
    if (token.empty()) {
      token = parser.argName();
    }
    // Set all forces to be the chosen force.
    gflow->forceMaster->setInteraction(ParseConstructor::getInteraction(parser.getNode(), variables, token, gflow));
  }
    // The interaction grid is specified.
  else {
    std::map<string, shared_ptr<Interaction> > interactions;
    string token;
    int t1, t2;
    parser.begin();
    // Go through each entry in the force table (each head is an entry).
    do {
      // CASE: [ : type1, type2, interaction-type ]
      if (parser.args_size() == 3) {
        // This allows the force string to be a variable.
        token = parser.get_variable_string(parser.argName(2));
        if (token.empty()) {
          token = parser.argName(2);
        }
        // Get the particle types
        t1 = parser.arg_cast<int>(0);
        t2 = parser.arg_cast<int>(1);
        // Make sure types are fine.
        if (t1 < 0 || t2 < 0 || NTypes <= t1 || NTypes <= t2) {
          throw BadStructure("Illegal particle type in force grid.");
        }
        // If the interaction has not occured yet, create one.
        if (interactions.find(token) == interactions.end()) {
          interactions.insert(make_pair(token,
                                        ParseConstructor::getInteraction(parser.getNode(), variables, token, gflow)));
        }
        // Set the interparticle interaction
        gflow->forceMaster->setInteraction(t1, t2, interactions.find(token)->second);
      }
        // NOT A VALID CASE
      else {
        throw BadStructure("Force grid: unrecognized structure.");
      }
    } while (parser.next());
    // We have scanned all the entries.
  }
  // Return to the original level.
  parser.up();

  // --- Global Particle Template. Defines "types" of particles, e.g. radius distribution, density/mass, etc.
  if (parser.begin("Template")) {
    do {
      ParseConstructor::parse_particle_template(parser.getNode(), variables, global_templates, gflow);
    } while (parser.next());
  }

  // --- Fill an area with particles
  if (parser.begin("Fill")) {
    FillAreaCreator fa(global_templates);
    PolymerCreator pc(global_templates);
    WallCreator wc(global_templates);
    CircleCreator cc(global_templates);
    do {
      if (parser.argName() == "Area") {
        fa.createArea(parser.getNode(), gflow, variables, particle_fixers);
      }
      else if (parser.argName() == "Polymer") {
        pc.createArea(parser.getNode(), gflow, variables, particle_fixers);
      }
      else if (parser.argName() == "Wall") {
        wc.createArea(parser.getNode(), gflow, variables, particle_fixers);
      }
      else if (parser.argName() == "Circle") {
        cc.createArea(parser.getNode(), gflow, variables, particle_fixers);
      }
      else {
        throw BadStructure("Unrecognized fill option: [" + parser.argName() + "].");
      }
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

  // Initialize handler
  gflow->handler->initialize();

  // --- Hard sphere relax the simulation. We do this by default.
  RealType relax_time = 0.5;
  parser.firstArg("HSRelax", relax_time);
  if (relax_time > 0) {
    hs_relax(gflow, relax_time);
  }

  // --- Natural relax of the simulation
  relax_time = 0.;
  parser.firstArg("Relax", relax_time);
  if (relax_time > 0) {
    relax(gflow, relax_time, parser.getNode("Relax"));
  }

  // --- Reconcile overlaps
  if (parser.begin("Reconcile")) {
    do {
      if (parser.argName() == "Remove") {
        // Default remove overlap cutoff.
        RealType rem = 2.;
        // Remove overlap cutoff
        parser.val(rem);
        // Do the removal
        gflow->removeOverlapping(rem);
      }
      else {
        throw BadStructure("Unrecognized remove option, [" + parser.argName() + "].");
      }
    } while (parser.next());
  }

  real skin = 0;
  if (parser.firstArg("SkinDepth", skin)) {
    gflow->handler->setSkinDepth(skin);
  }

  // Reconstruct handler.
  gflow->handler->construct();

  // Fix velocities
  fix_particle_velocities(gflow->simData);
}

inline Integrator *FileParseCreator::choose_integrator(HeadNode *head) const {
  // Create parser, with variables.
  TreeParser parser(head, variables);
  // Declare valid options.
  parser.addHeadingOptional("Number");
  parser.addHeadingOptional("MaxDT");
  parser.addHeadingOptional("MinDT");
  parser.addHeadingOptional("Adjust");
  parser.addHeadingOptional("UseV");
  parser.addHeadingOptional("UseA");

  int ivalue = 0;
  RealType rvalue = 0.;
  bool bvalue = false;

  string token = parser.argName();

  Integrator *integrator = nullptr;
  if (token == "VelocityVerlet") {
    integrator = choose_velocity_verlet(gflow, sim_dimensions);
  }
  else if (token == "OverdampedIntegrator") {
    integrator = choose_overdamped_integrator(gflow, sim_dimensions);
  }
  else if (token == "OverdampedLangevin") {
    integrator = new OverdampedLangevinIntegrator(gflow);
  }
  else if (token == "LangevinIntegrator") {
    integrator = new LangevinIntegrator(gflow);
  }
  else if (token == "NoseHooverVelocityVerlet") {
    integrator = choose_nose_hoover_velocity_verlet(gflow, sim_dimensions);
  }
  else {
    throw UnexpectedOption("Integrator choice was [" + token + "].");
  }
  // Temperature and viscosity for LangevinType integrators
  auto lti = dynamic_cast<LangevinTypeIntegrator *>(integrator);

  if (lti) {
    parser.addHeadingOptional("Temperature");
    parser.addHeadingOptional("Viscosity");
    parser.check();
    // Look for temperature or viscosity
    if (parser.firstArg("Temperature", rvalue)) {
      lti->setTemperature(rvalue);
    }
    if (parser.firstArg("Viscosity", rvalue)) {
      lti->setViscosity(rvalue);
    }
  }
  else {
    parser.check();
  }

  // Options for any integrator
  if (parser.firstArg("Number", ivalue)) {
    integrator->setStepDelay(ivalue);
  }
  if (parser.firstArg("MaxDT", rvalue) && rvalue > 0) {
    integrator->setMaxDT(rvalue);
  }
  if (parser.firstArg("MinDT", rvalue) && rvalue > 0) {
    integrator->setMinDT(rvalue);
  }
  if (parser.firstArg("Adjust", bvalue)) {
    integrator->setAdjustDT(bvalue);
  }
  if (parser.firstArg("UseV", bvalue)) {
    integrator->setUseV(bvalue);
  }
  if (parser.firstArg("UseA", bvalue)) {
    integrator->setUseA(bvalue);
  }

  // Return
  return integrator;
}

inline shared_ptr<Interaction> FileParseCreator::choose_interaction(HeadNode *head) const {
  string token = head->params[2]->partA;
  return InteractionChoice::choose(gflow, token, sim_dimensions);
}

inline BCFlag FileParseCreator::choose_bc(const string &token) const {
  if (token == "Wrap") {
    return BCFlag::WRAP;
  }
  else if (token == "Reflect") {
    return BCFlag::REFL;
  }
  else if (token == "Repulse") {
    return BCFlag::REPL;
  }
  else if (token == "Open") {
    return BCFlag::OPEN;
  }
  else {
    throw UnexpectedOption("Boundary condition choice was [" + token + "].");
  }
}

inline void FileParseCreator::add_modifier(HeadNode *head) const {
  TreeParser parser(head, variables);
  string token = parser.argName();

  if (token == "WindTunnel") {
    auto wind_tunnel = make_shared<WindTunnel>(gflow);
    wind_tunnel->parse_construct(head, variables);
    gflow->addModifier(wind_tunnel);
  }
  else if (token == "StreamTunnel") {
    auto stream_tunnel = make_shared<StreamTunnel>(gflow);
    stream_tunnel->parse_construct(head, variables);
    gflow->addModifier(stream_tunnel);
  }
  else {
    throw UnexpectedOption("Modifier choice was [" + token + "].");
  }
}

inline void FileParseCreator::createParticle(HeadNode *head) const {
  // Check if head is good
  if (head == nullptr) {
    return;
  }

  TreeParser parser(head, variables);
  // Declare valid options.
  parser.addHeadingNecessary("Position", "Particle needs position data.");
  parser.addHeadingNecessary("Velocity", "Particle needs velocity data.");
  parser.addHeadingNecessary("Sigma", "Particle needs radius data.");
  parser.addHeadingOptional("Type");
  parser.addHeadingOptional("Mass");
  parser.addHeadingOptional("Modifier");
  // Check headings for validity.
  parser.check();

  // --- Look for options
  Vec X(sim_dimensions), V(sim_dimensions);
  RealType sigma, im;
  int type = 0;

  parser.firstArgVec("Position", X);
  parser.firstArgVec("Velocity", V);
  parser.firstArg("Sigma", sigma);
  parser.firstArg("Type", type);

  // Get mass
  if (parser.focus("Mass")) {
    RealType d = 1;
    string token = parser.argName();
    if (token == "inf") {
      im = 0;
    }
    else if (token == "Density") {
      if (parser.val(d)) {
        im = 1. / (d * sphere_volume(sigma, sim_dimensions));
      }
      else {
        throw BadStructure("Need to specify the density.");
      }
    }
    else { // It must be the mass
      if (parser.arg(d)) {
        im = 1. / d;
      }
      else {
        throw BadStructure("Need to specify the mass.");
      }
    }
  }
  else {
    im = 1. / sphere_volume(sigma, sim_dimensions);
  } // Default density is 1
  // Return to top level.
  parser.up();

  // If there is a modifier, set that up.
  if (parser.begin("Modifier")) {
    // Next global id - this will be the added particle's global id.
    int g_id = gflow->simData->getNextGlobalID();
    // Go through all modifiers.
    do {
      string token = parser.argName();
      // Check token
      if (token == "CV") {
        gflow->addModifier(make_shared<ConstantVelocity>(gflow, g_id, V.data));
      }
      else if (token == "CV-D") {
        RealType D = 0;
        if (parser.arg(D)) {
          gflow->addModifier(make_shared<ConstantVelocityDistance>(gflow, g_id, V.data, D));
        }
        else {
          throw BadStructure("Distance needs to be specified for CV-D modifier.");
        }
      }
      else {
        throw BadStructure("Unrecognized modifer option, [" + token + "].");
      }
    } while (parser.next());
  }

  // Add the particle to the system
  gflow->simData->addParticle(X.data, V.data, sigma, im, type);
}

inline void FileParseCreator::makeRandomForces() {
  // Assign random interactions, either HardSphere, LennardJones, or Coulomb (for now), and with equal probability (for now)
  auto hs = InteractionChoice::choose(gflow, HardSphereToken, sim_dimensions);
  auto lj = InteractionChoice::choose(gflow, LennardJonesToken, sim_dimensions);
  auto cl = InteractionChoice::choose(gflow, CoulombToken, sim_dimensions);
  // Assign random (but symmetric) interactions
  for (int i = 0; i < NTypes; ++i) {
    // Self interaction
    float choice = drand48();
    if (choice < 0.333) {
      gflow->forceMaster->setInteraction(i, i, hs);
    }
    else if (choice < 0.666) {
      gflow->forceMaster->setInteraction(i, i, lj);
    }
    else {
      gflow->forceMaster->setInteraction(i, i, cl);
    }
    // Interactions with other objects.
    for (int j = i + 1; j < NTypes; ++j) {
      choice = drand48();
      if (choice < 0.333) {
        gflow->forceMaster->setInteraction(i, j, hs);
      }
      else if (choice < 0.666) {
        gflow->forceMaster->setInteraction(i, j, lj);
      }
      else {
        gflow->forceMaster->setInteraction(i, j, cl);
      }
    }
  }
}

inline string FileParseCreator::copyFile() const {
  string s;
  char c;
  ifstream fin(configFile);
  if (fin.fail()) {
    return "";
  }
  fin.get(c);
  while (!fin.eof()) {
    s.push_back(c);
    fin.get(c);
  }
  return s;
}
