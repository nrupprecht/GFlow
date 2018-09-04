#include "fileparsecreator.hpp"
// Other files
#include "../utility/printingutility.hpp"
#include "../allmodifiers.hpp"

namespace GFlowSimulation {

  FileParseCreator::FileParseCreator(int argc, char **argv) : Creator(argc, argv), configFile(""), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    // Seed generators here
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  FileParseCreator::FileParseCreator(ArgParse *p) : Creator(p), configFile(""), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    // Seed generators here
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  FileParseCreator::FileParseCreator(ArgParse *p, string f) : Creator(p), configFile(f), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    // Seed generators here
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  GFlow* FileParseCreator::createSimulation() {
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

    // Sort and collect options
    build_message += "Collecting options... ";
    std::multimap<string, HeadNode*> options;
    for (const auto &op : root->subHeads) {
      options.insert(std::pair<string, HeadNode*> (op->heading, op));
    }
    build_message += "Done.\n";

    // Create the scenario from the options
    if (gflow) delete gflow;
    gflow = new GFlow;
    build_message += "Starting simulation setup... ";
    try {
      createFromOptions(gflow, options);
    }
    catch (BadStructure bs) {
      cout << build_message << endl;
      cout << "Caught Bad Structure error: " << bs.message << endl;
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

    // Clean up and return
    //delete [] root;
    return gflow;
  }

  inline void FileParseCreator::createFromOptions(GFlow *gflow, std::multimap<string, HeadNode*>& options) {
    // Container to collect options in
    std::vector<HeadNode*> container;
    // A head node we can use
    HeadNode *head;

    // --- Look for dimensions
    getAllMatches("Dimensions", container, options);
    if (container.size()>1) build_message +=  "We only need the number of dimensions specified once. Ignoring all but the first.\n";
    if (container.size()==1) {
      head = container[0];
      if (head->params.size()!=1) 
        throw BadStructure("Dimensions should be a single argument, we found "+toStr(head->params.size()));
      if (convert<int>(head->params[0]->partA)!=DIMENSIONS) throw BadDimension();
    }
    
    // --- Look for bounds
    getAllMatches("Bounds", container, options);
    if (container.size()==0) throw BadStructure("We need bounds information.");
    if (container.size()>1) build_message +=  "We only need the number of dimensions specified once. Ignoring all but the first.\n";
    head = container[0];
    if (head->subHeads.size()!=DIMENSIONS) 
      throw BadStructure("For bounds, we need "+toStr(DIMENSIONS)+" conditions, we found "+toStr(head->subHeads.size()));
    for (int d=0; d<head->subHeads.size(); ++d) {
      if (head->subHeads[d]->params.size()!=2) 
        throw BadStructure("Bounds need a min and a max, we found "+toStr(head->subHeads[d]->params.size()+" parameters."));
      bounds.min[d] = convert<float>( head->subHeads[d]->params[0]->partA );
      bounds.max[d] = convert<float>( head->subHeads[d]->params[1]->partA );
    }
    gflow->setBounds(bounds);

    // --- Look for integrator
    getAllMatches("Integrator", container, options);
    if (container.size()>1) build_message +=  "We only need the integrator type specified once. Ignoring all but the first.\n";
    // Default integrator is velocity verlet
    if (container.empty()) gflow->integrator = new VelocityVerlet(gflow);
    else {
      head = container[0];
      // We expect a single option: the name of the type of integrator to use
      if (head->params.size()!=1) throw BadStructure("In Integrator we found more than one word.");
      // Set new integrator
      gflow->integrator = choose_integrator(head);
    }
    
    // --- Look for number of particle types
    getAllMatches(Types_Token, container, options);
    if (container.size()>1) build_message += "We only need the number of types specified once. Ignoring all but the first.\n";
    // Default is 1 type
    if (container.empty()) NTypes = 1;
    else {
      head = container[0];
      if (head->params.size()!=1) 
        throw BadStructure("In NTypes we found more than one number.");
      // Set particle types
      NTypes = convert<int>(head->params[0]->partA);
      gflow->forceMaster->setNTypes(NTypes);
    }

    // --- Look for interactions
    getAllMatches(Interactions_Token, container, options);
    if (container.size()>1) build_message += "We only need the interactions specified once. Ignoring all but the first.\n";
    head = container[0];
    // Collect all the interactions we need
    std::map<string, Interaction*> interactions;
    // The body specifies the force grid
    for (auto fg : head->subHeads) {
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

    // --- Look for boundary conditions
    getAllMatches("Boundary-conditions", container, options);
    if (container.size()>1) build_message += "We only need the boundary conditions specified once. Ignoring all but the first.\n";
    if (container.size()>0) {
      head = container[0];
      // If a parameter is given, it is the BC for all sides
      if (head && !head->params.empty()) gflow->setAllBCs(choose_bc(head->params[0]->partA));
      else {
        for (auto bc : head->subHeads) {
          if (!head->heading.empty() || head->params.size()==1) { // A specific dimension must be choosen, the (one) parameter is the flag
            if (head->params.empty()) throw BadStructure("Expected a parameter, found none.");
            gflow->setBC(convert<int>(head->heading), choose_bc(head->params[0]->partA));
          }
          else throw BadStructure("We need one parameter to be the boundary condition flag, we found "+toStr(head->params.size()));
        }
      }
    }

    // --- Global Particle Template. Defines "types" of particles, e.g. radius distribution, density/mass, etc.
    getAllMatches("Template", container, options);
    for (auto h : container) {
      getParticleTemplate(h, global_templates);
    }    

    getAllMatches("Fill-area", container, options);
    for (auto h : container) {
      fillArea(h);
    }

    getAllMatches("Creation", container, options);
    if (container.size()>1) build_message += "We only want one specification of particle creation. Ignoring all but the first";
    if (container.size()>0) {
      head = container[0];
      // The subheads have no bodies, and are of the form: [type i] : [type j], [rate].
      // For now, we only allow reproduction to create the same type.
      std::map<int, RealType> birth;
      for (auto s : head->subHeads) {
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

    getAllMatches("Destruction", container, options);
    if (container.size()>1) build_message += "We only want one specification of particle destruction. Ignoring all but the first";
    if (container.size()>0) {
      head = container[0];
      // The subheads have no bodies, and are of the form: [type i] : [rate].
      // For now, we only allow reproduction to create the same type.
      std::map<int, RealType> death;
      for (auto s : head->subHeads) {
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

    // Initialize domain
    gflow->domain->initialize();
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
    else if (token=="LennardJones") return new LennardJones(gflow);
    else if (token=="Coulomb")      return new CoulumbForce(gflow);
    else if (token=="Consumption")  return new Consumption(gflow);
    else if (token=="None")         return nullptr;
    else throw UnexpectedOption();
  }

  inline BCFlag FileParseCreator::choose_bc(string& token) const {
    if      (token=="Wrap") return BCFlag::WRAP;
    else if (token=="Refl") return BCFlag::REFL;
    else if (token=="Repl") return BCFlag::REPL;
    else throw UnexpectedOption();
  }

  inline void FileParseCreator::fillArea(HeadNode *head) const {
    // Check if head is good
    if (head==nullptr) return;
    // Particle Templates
    std::map<string, ParticleTemplate> particle_templates;

    // --- Sort options
    std::multimap<string, HeadNode*> options;
    for (auto h : head->subHeads)
      options.insert(std::pair<string, HeadNode*>(h->heading, h));
    // Container to collect options in
    std::vector<HeadNode*> container;

    // --- Look for options

    // --- Look for bounds
    getAllMatches("Bounds", container, options);
    // We must have bounds
    Bounds bnds; 
    // We only care about the first bounds. We can only define the bounds once.
    if (container.empty())
      throw BadStructure("We need bounds!");
    else {
      // We want there to only be one set of bounds
      if (container.size()>1) build_message += "Multiple bounds found! Ignoring all but the first instance.\n";
      // Get the bounds from the first instance
      HeadNode *h = container[0];
      // Check if we should use the full simulation bounds
      if (h->subHeads.empty() && h->params.size()>0 && h->params[0]->partA=="Full") {
        bnds = bounds;
      }
      else {
        if (h->subHeads.size()!=DIMENSIONS) 
          throw BadStructure("Expected "+toStr(DIMENSIONS)+" arguments, found "+toStr(h->subHeads.size()));
        // Set bounds
        for (int d=0; d<h->subHeads.size(); ++d) {
          // Check for well formed options.
          if (h->subHeads[d]->params.size()!=2 || !h->subHeads[d]->subHeads.empty()) 
            throw BadStructure("Bounds need a min and a max, we found "+toStr(h->subHeads[d]->params.size())+" parameters.");
          // Extract the bounds.
          bnds.min[d] = convert<float>( h->subHeads[d]->params[0]->partA );
          bnds.max[d] = convert<float>( h->subHeads[d]->params[1]->partA );
        }
      }
    }

    // --- Local Particle Template. Defines "types" of particles, e.g. radius distribution, density/mass, etc.
    getAllMatches("Template", container, options);
    for (auto h : container) {
      getParticleTemplate(h, particle_templates);
    }

    // --- Number. How to choose which particles to fill the space with.
    getAllMatches("Number", container, options);
    if (container.empty()) 
      throw BadStructure("We need number information!");
    // Create a structure for recording probabilities or numbers
    std::map<string, double> particle_template_numbers;
    bool useNumber = false, usePhi = false, singleType = false;
    int number(0); 
    double phi(0);
    // Find options
    if (container.size()>1) build_message += "Multiple number directives found. Ignoring all but the first instance.\n";
    HeadNode *hd = container[0];
    // No body - we must have something in one of two forms:
    // Number: #
    // Number: Phi=#
    if (hd->subHeads.size()==0) { 
      if (particle_templates.size()>1) 
        throw BadStructure("More than one type of particle has been defined, but how probable they are is ill defined.");
      if (hd->params.empty())
        throw BadStructure("Expected parameters in number directive, since there is no body.");
      // Single type scenario
      singleType = true;
      // Phi or number?
      if (hd->params[0]->partB.empty()) {
        useNumber = true;
        number = convert<int>(hd->params[0]->partA);
      }
      else if (hd->params[0]->partA=="Phi") {
        usePhi = true;
        phi = convert<double>(hd->params[0]->partB);
      }
      else throw BadStructure("Expect either a number or 'Phi=#'");
    }
    // Mutiple particle templates are defined. There must be a body
    else {
      // Either phi or number
      if (!hd->params.empty() && hd->params[0]->partA=="Phi") usePhi = true;
      else useNumber = true;
      // Look in the body for the probabilities or numbers
      for (auto sh : hd->subHeads) {
        if (sh->params.size()!=1)
          throw BadStructure("Expect two options in template - number definition. Found "+toStr(sh->params.size())+".");
        // Insert numbers
        particle_template_numbers.insert(
          std::pair<string, double>(sh->heading, convert<double>(sh->params[0]->partA))
        );
      }
    }

    // Select a velocity
    auto select_velocity = [&] (RealType *V, RealType *X, RealType sigma, RealType im, int type) -> void {
      RealType vsgma = 0.25;
      // Velocity based on KE
      double ke = fabs(vsgma*normal_dist(generator));
      double velocity = sqrt(2*im*ke/127.324);
      // Random normal vector
      RealType normal[DIMENSIONS];
      randomNormalVec(normal);
      // Set the velocity
      scalarMultVec(velocity, normal, V);
    };    

    // --- Check that we have defined a good area
    for (int d=0; d<DIMENSIONS; ++d)
      if (bnds.wd(d)==0) 
	 throw BadStructure("We need valid bounds. The width of dimension "+toStr(d)+" was 0.");
    if (!usePhi && number<=0 && singleType)
      throw BadStructure("If using a single type, we need a nonzero number of particles.");

    // --- We have found all the options. Fill the area.
    GFlow filler;
    filler.setBounds(bnds);
    filler.setAllBCs(BCFlag::WRAP);
    filler.forceMaster = gflow->forceMaster; // Make sure the particles treat each other in the same way
    // Get the simdata
    SimData *simData = filler.simData;
    // --- Fill with particles
    RealType X[DIMENSIONS], V[DIMENSIONS], sigma(0.), im(0.);
    int type(0);
    zeroVec(V);

    // If we are filling to a specified packing fraction
    if (usePhi) {
      // One type of particle
      if (phi>0) {
        // We expect there to only be one particle template
        if (particle_templates.size()!=1)
          throw BadStructure("We expect there to be 1 particle template, there are "+toStr(particle_templates.size()));
        // Get the particle template
        ParticleTemplate &particle_creator = particle_templates.begin()->second;
        int i(0);
        RealType vol = 0, Vol = bnds.vol();
        while (vol/Vol < phi) {
          // Select a position for the particle (random uniform)
          for (int d=0; d<DIMENSIONS; ++d)
            X[d] = drand48()*bnds.wd(d) + bnds.min[d];
          // Select other characteristics
          particle_creator.createParticle(X, sigma, im, type, i);
          // Add the particle
          simData->addParticle(X, V, sigma, im, type);
          // Increment volume and counter
          vol += sphere_volume(sigma);
          ++i;
        }
      }
      // Multiple types of particles, each must have different target phi's
      else {
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
        RealType vol = 0, Vol = bnds.vol();
        while (vol/Vol < phi) {
          // Select a position for the particle (random uniform)
          for (int d=0; d<DIMENSIONS; ++d)
            X[d] = drand48()*bnds.wd(d) + bnds.min[d];
          // Choose a type of particle to create
          int pt = choice(global_generator);
          ParticleTemplate &particle_creator = template_vector.at(pt);
          // Select other characteristics
          particle_creator.createParticle(X, sigma, im, type, i);
          // Add the particle
          simData->addParticle(X, V, sigma, im, type);
          // Increment volume and counter
          vol += sphere_volume(sigma);
          ++i;
        }
      }
    }
    // If we are filling to a specified number
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
          for (int d=0; d<DIMENSIONS; ++d)
            X[d] = drand48()*bnds.wd(d) + bnds.min[d];
          // Select other characteristics
          particle_creator.createParticle(X, sigma, im, type, i);
          // Add the particle
          simData->addParticle(X, V, sigma, im, type);
        }
      }
    }

    // Print status
    build_message += "From Fill Area: Done with initial particle assigmnemt.\n"; 

    // Initialize domain
    filler.domain->initialize();
    filler.integrator = new VelocityVerlet(&filler);

    // --- Relax the simulation
    hs_relax(&filler, 0.1); // 1) To make sure particles don't stop on top of one another
    relax(&filler, 0.15);

    // --- Fill gflow with the particles
    for (int i=0; i<simData->number; ++i) {
      // Extract the particle properties
      copyVec(simData->X(i), X);
      int type = simData->Type(i);
      RealType sigma = simData->Sg(i);
      RealType im = simData->Im(i);
      if (type!=-1) {
        // Select the velocity for the final particle
        select_velocity(V, X, sigma, im, type);
        if (im==0) zeroVec(V); //** FOR NOW
        gflow->simData->addParticle(X, V, sigma, im, type);
      }
    }
      
    // So we don't delete the force master when filler cleans up
    filler.forceMaster = nullptr; 
  }

  inline void FileParseCreator::getAllMatches(string heading, vector<HeadNode*>& container, std::multimap<string, HeadNode*>& options) const {
    // Clear container in case there is stuff left over in it.
    container.clear();
    // Look for options
    bool good = true;
    for (auto it=options.find(heading); it!=options.end() && good; ++it) {
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
    p_template.mass_engine = getRandomEngine(container[0], type);
    p_template.mass_string = type;

    // --- Look for Type option
    getAllMatches("Type", container, options);
    if (container.empty())
      throw BadStructure("We need some mass information!");
    if (container.size()>1)
      build_message += "We only need one mass information block. Ignoring all but the first instance.\n";
    p_template.type_engine = getRandomEngine(container[0], type);
    p_template.type_string = type;

    // Add to particle templates
    particle_templates.insert(std::pair<string, ParticleTemplate>(head->params[0]->partA, p_template));
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
    else if (h->subHeads.empty()) {
      string token = h->params[0]->partA;
      // We expect either a number in partA, or a string (e.g. "Inf").
      if (isdigit(token.at(0))) {
        return new DeterministicEngine(convert<double>(token));
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
      if (lwr->heading!="Lower") std::swap(lwr, upr);
      // Check structure
      if (lwr->heading!="Lower" || upr->heading!="Upper") 
        throw BadStructure("Headings incorrect for uniform distribution.");
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

}
