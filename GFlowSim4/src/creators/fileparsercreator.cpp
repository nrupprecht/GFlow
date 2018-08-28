#include "fileparsecreator.hpp"
// Other files
#include "../utility/printingutility.hpp"

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
    cout << "Starting file parse... ";
    FileParse parser;
    HeadNode *root = parser.parseFile(configFile); // A file parse class does the parsing
    cout << "Done.\n";

    // Get the message from the parser.
    parse_message = parser.getMessage();

    // Sort and collect options
    cout << "Collecting options... ";
    std::multimap<string, HeadNode*> options;
    for (const auto &op : root->subHeads) {
      options.insert(std::pair<string, HeadNode*> (op->heading, op));
    }
    cout << "Done collecting options.\n";

    // Create the scenario from the options
    gflow = new GFlow;
    cout << "Starting simulation setup... ";
    try {
      createFromOptions(gflow, options);
    }
    catch (BadStructure bs) {
      cout << bs.message << endl;
      throw;
    }
    cout << "Done.\n";

    // Clean up and return
    //delete [] root;
    return gflow;
  }

  inline void FileParseCreator::createFromOptions(GFlow *gflow, std::multimap<string, HeadNode*>& options) const {
    string token;
    bool good = true;
    std::multimap<string, HeadNode*>::iterator it;

    // Look for dimensions
    good = true;
    token = "Dimensions";
    it = options.find(token);
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // For now, we can only check whether the dimensions are consistent, DIMENSIONS is a global constant.
        if (h->params.size()!=1) throw BadStructure("Dimensions should be a single argument, we found "+toStr(h->params.size()));
        if (convert<int>(h->params[0]->partA)!=DIMENSIONS) throw BadDimension();
      }
      else good = false;
    }

    // Look for bounds
    token = "Bounds";
    good = true;
    it = options.find(token);
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        Bounds bounds;
        // The body of bounds contains the actual bounds
        if (h->subHeads.size()!=DIMENSIONS) 
          throw BadStructure("For bounds, we need "+toStr(DIMENSIONS)+" conditions, we found "+toStr(h->subHeads.size()));
        // Set bounds
        for (int d=0; d<h->subHeads.size(); ++d) {
          if (h->subHeads[d]->params.size()!=2) 
            throw BadStructure("Bounds need a min and a max, we found "+toStr(h->subHeads[d]->params.size()+" parameters."));
          bounds.min[d] = convert<float>( h->subHeads[d]->params[0]->partA );
          bounds.max[d] = convert<float>( h->subHeads[d]->params[1]->partA );
        }
        gflow->setBounds(bounds);
      }
      else good = false;
    }

    // Look for Integrator
    token = "Integrator";
    good = true;
    it = options.find(token);
    if (it==options.end()) { // No integrator was specified. Use a velocity verlet integrator.
      gflow->integrator = new VelocityVerlet(gflow);
    }
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // We expect a single option: the name of the type of integrator to use
        if (h->params.size()!=1) throw BadStructure("In Integrator we found more than one word.");
        // Get rid of old integrator if necessary
        if (gflow->integrator) delete gflow->integrator;
        // Set new integrator
        gflow->integrator = choose_integrator(h->params[0]->partA);
      }
      else good = false;
    }

    // Look for number of particle types
    token = "NTypes";
    good = true;
    it = options.find(token);
    if (it==options.end()) { // No integrator was specified. Use a velocity verlet integrator.
      gflow->integrator = new VelocityVerlet(gflow);
    }
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // We expect a single option: the number of particle types
        if (h->params.size()!=1) 
          throw BadStructure("In NTypes we found more than one number.");
        // Set particle types
        gflow->forceMaster->setNTypes(convert<int>(h->params[0]->partA));
      }
      else good = false;
    }

    // Look for a force grid
    token = "Force-grid";
    good = true;
    it = options.find(token);
    if (it==options.end()) { // No integrator was specified. Use a velocity verlet integrator.
      gflow->integrator = new VelocityVerlet(gflow);
    }
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // Collect all the interactions we need
        std::map<string, Interaction*> interactions;
        // The body specifies the force grid
        for (auto fg : h->subHeads) {
          // Look for type1, type2, interaction-type
          if (fg->params.size()!=3) 
            throw BadStructure("Force grid needs three parameters per line to specify interaction, found "+toStr(fg->params.size())+".");
          // What type of interaction do we need
          string i_token = fg->params[2]->partA;
          if (interactions.find(i_token)==interactions.end()) { // New interaction type
            interactions.insert(std::pair<string, Interaction*>(i_token, choose_interaction(i_token)));
          }
        }
        // Now assign types to the interactions we have found. We have already checked for structure, no need to check again.
        for (auto fg : h->subHeads) {
          int t1 = convert<int>(fg->params[0]->partA), t2 = convert<int>(fg->params[1]->partA);
          string i_token = fg->params[2]->partA;
          gflow->forceMaster->setInteraction(t1, t2, interactions.find(i_token)->second);
        }
      }
      else good = false;
    }

    // Look for boundary conditions
    token = "Boundary-conditions";
    good = true;
    it = options.find(token);
    // Default option is all wrapped boundaries
    gflow->setAllBCs(BCFlag::WRAP);
    // Look for options
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // If a parameter is given, it is the BC for all sides
        if (h && !h->params.empty()) gflow->setAllBCs(choose_bc(h->params[0]->partA));
        else {
          for (auto bc : h->subHeads) {
            if (!h->heading.empty() || h->params.size()==1) { // A specific dimension must be choosen, the (one) parameter is the flag
              gflow->setBC(convert<int>(h->heading), choose_bc(h->params[0]->partA));
            }
            else throw BadStructure("We need one parameter to be the boundary condition flag, we found "+toStr(h->params.size()));
          }
        }
      }
      else good = false;
    }

    // Look for fill areas
    token = "Fill-area";
    good = true;
    it = options.find(token);
    // Look for options
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // This is complicated enough that we give it it's own function. Fill area has its own options.
        if (h) fillArea(h);
      }
      else good = false;
    }

    // Initialize domain
    gflow->domain->initialize();
  }

  inline Integrator* FileParseCreator::choose_integrator(string& token) const {
    if (token=="VelocityVerlet")            return new VelocityVerlet(gflow);
    else if (token=="OverdampedIntegrator") return new OverdampedIntegrator(gflow);
    else throw UnexpectedOption();
  }

  inline Interaction* FileParseCreator::choose_interaction(string& token) const {
    if (token=="HardSphere")        return new HardSphere(gflow);
    else if (token=="LennardJones") return new LennardJones(gflow);
    else throw UnexpectedOption();
  }

  inline BCFlag FileParseCreator::choose_bc(string& token) const {
    if (token=="Wrap")      return BCFlag::WRAP;
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
    Bounds bounds; 
    // We only care about the first bounds. We can only define the bounds once.
    if (container.empty())
      throw BadStructure("We need bounds!");
    else {
      // We want there to only be one set of bounds
      if (container.size()>1) cout << "Multiple bounds found! Ignoring all but the first instance.\n";
      // Get the bounds from the first instance
      HeadNode *h = container[0];
      if (h->subHeads.size()!=DIMENSIONS) 
        throw BadStructure("Expected "+toStr(DIMENSIONS)+" arguments, found "+toStr(h->subHeads.size()));
      // Set bounds
      for (int d=0; d<h->subHeads.size(); ++d) {
        // Check for well formed options.
        if (h->subHeads[d]->params.size()!=2 || !h->subHeads[d]->subHeads.empty()) 
          throw BadStructure("Bounds need a min and a max, we found "+toStr(h->subHeads[d]->params.size())+" parameters.");
        // Extract the bounds.
        bounds.min[d] = convert<float>( h->subHeads[d]->params[0]->partA );
        bounds.max[d] = convert<float>( h->subHeads[d]->params[1]->partA );
      }
    }

    // --- Particle Template. Defines "types" of particles, e.g. radius distribution, density/mass, etc.
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
    bool useNumber = false, usePhi = false;
    int number; 
    double phi;
    // Find options
    if (container.size()>1) cout << "Multiple number directives found. Ignoring all but the first instance.\n";
    HeadNode *hd = container[0];
    if (hd->subHeads.size()==0) { // No body
      if (particle_templates.size()>1) 
        throw BadStructure("More than one type of particle has been defined, but how probable they are is ill defined.");
      if (hd->params.empty())
        throw BadStructure("Expected parameters in number directive, since there is no body.");
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
      if (!hd->params.empty() && hd->params[0]->partA=="Phi") {
        usePhi = true;
        phi = convert<double>(hd->params[0]->partB);
      }
      else useNumber = true;
      // Look in the body for the probabilities or frequencies
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
      if (bounds.wd(d)==0) return;
    if (number<=0 && !usePhi) return;

    // --- We have found all the options. Fill the area.
    GFlow filler;
    filler.setBounds(bounds);
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
      RealType vol = 0, Vol = bounds.vol();
      while (vol/Vol < phi) {
        // Select a position for the particle (random uniform)
        for (int d=0; d<DIMENSIONS; ++d)
          X[d] = drand48()*bounds.wd(d) + bounds.min[d];
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
            X[d] = drand48()*bounds.wd(d) + bounds.min[d];
          // Select other characteristics
          particle_creator.createParticle(X, sigma, im, type, i);
          // Add the particle
          simData->addParticle(X, V, sigma, im, type);
        }
      }
    }

    // Print status
    cout << "Done with initial particle assigmnemt.\n"; 

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
    if (container.size()>1) {
      cout << "We only need one radius information block. Ignoring all but the first instance.\n";
    }
    // Get the random engine
    p_template.radius_engine = getRandomEngine(container[0], type);
    
      
    // --- Look for Mass option
    getAllMatches("Mass", container, options);
    if (container.empty())
      throw BadStructure("We need some mass information!");
    if (container.size()>1) {
      cout << "We only need one mass information block. Ignoring all but the first instance.\n";
    }
    p_template.mass_engine = getRandomEngine(container[0], type);

    // --- Look for Type option
    getAllMatches("Type", container, options);
    if (container.empty())
      throw BadStructure("We need some mass information!");
    if (container.size()>1) {
      cout << "We only need one mass information block. Ignoring all but the first instance.\n";
    }
    p_template.type_engine = getRandomEngine(container[0], type);

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
      // We expect a number in partA
      return new DeterministicEngine(convert<double>(h->params[0]->partA));
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
        probabilities.push_back(convert<double>(sh->params[1]->partB));
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