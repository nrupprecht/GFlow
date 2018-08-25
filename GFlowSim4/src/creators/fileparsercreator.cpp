#include "fileparsecreator.hpp"
// Other files
#include "../utility/printingutility.hpp"

namespace GFlowSimulation {

  FileParseCreator::FileParseCreator(int argc, char **argv) : Creator(argc, argv), configFile(""), root(nullptr), currentHead(nullptr), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    // Seed generators here
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  FileParseCreator::FileParseCreator(ArgParse *p) : Creator(p), configFile(""), root(nullptr), currentHead(nullptr), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    // Seed generators here
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  FileParseCreator::FileParseCreator(ArgParse *p, string f) : Creator(p), configFile(f), root(nullptr), currentHead(nullptr), gflow(nullptr) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    seedGenerator(seed);
    // Seed generators here
    normal_dist = std::normal_distribution<RealType>(0., 1.);
  }

  GFlow* FileParseCreator::createSimulation() {
    // Create filestream
    std::ifstream fin(configFile);
    if (fin.fail()) {
      cout << "File parse creator failed to open file [" << configFile << "].\nExiting.\n";
      return nullptr;
    }

    // Set up head node
    root = new HeadNode;
    currentHead = root;

    // We treat the file as one giant body. Get that body.
    level = 0; 
    getBody(fin); // Parsing happens here

        // Sort and collect options
    std::multimap<string, HeadNode*> options;
    for (const auto &op : root->subHeads) {
      options.insert(std::pair<string, HeadNode*> (op->heading, op));
    }

    // Create the scenario from the options
    gflow = new GFlow;
    createFromOptions(gflow, options);
    return gflow;
  }

  inline void FileParseCreator::createFromOptions(GFlow *gflow, std::multimap<string, HeadNode*>& options) const {
    string token;
    bool good = true;
    std::multimap<string, HeadNode*>::iterator it;

    // Look for dimensions
    good = true;
    token = "Dimensions:";
    it = options.find(token);
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // For now, we can only check whether the dimensions are consistent, DIMENSIONS is a global constant.
        if (h->params.size()!=1) throw UnexpectedOption();
        if (convert<int>(h->params[0]->partA)!=DIMENSIONS) throw BadDimension();
      }
      else good = false;
    }

    // Look for bounds
    token = "Bounds:";
    good = true;
    it = options.find(token);
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        Bounds bounds;
        // The body of bounds contains the actual bounds
        if (h->subHeads.size()!=DIMENSIONS) throw BadStructure();
        // Set bounds
        for (int d=0; d<h->subHeads.size(); ++d) {
          if (h->subHeads[d]->params.size()!=2) throw BadStructure();
          bounds.min[d] = convert<float>( h->subHeads[d]->params[0]->partA );
          bounds.max[d] = convert<float>( h->subHeads[d]->params[1]->partA );
        }
        gflow->setBounds(bounds);
      }
      else good = false;
    }

    // Look for Integrator
    token = "Integrator:";
    good = true;
    it = options.find(token);
    if (it==options.end()) { // No integrator was specified. Use a velocity verlet integrator.
      gflow->integrator = new VelocityVerlet(gflow);
    }
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // We expect a single option: the name of the type of integrator to use
        if (h->params.size()!=1) throw BadStructure();
        // Get rid of old integrator if necessary
        if (gflow->integrator) delete gflow->integrator;
        // Set new integrator
        gflow->integrator = choose_integrator(h->params[0]->partA);
      }
      else good = false;
    }

    // Look for number of particle types
    token = "NTypes:";
    good = true;
    it = options.find(token);
    if (it==options.end()) { // No integrator was specified. Use a velocity verlet integrator.
      gflow->integrator = new VelocityVerlet(gflow);
    }
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // We expect a single option: the number of particle types
        if (h->params.size()!=1) throw BadStructure();
        // Set particle types
        gflow->forceMaster->setNTypes(convert<int>(h->params[0]->partA));
      }
      else good = false;
    }

    // Look for a force grid
    token = "Force-grid:";
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
          if (fg->params.size()!=3) throw BadStructure();
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
    token = "Boundary-conditions:";
    good = true;
    it = options.find(token);
    // Default option is all wrapped boundaries
    gflow->setAllBCs(BCFlag::WRAP);
    // Look for options
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // If a parameter is given, it is the BC for all sides
        if (!h->params.empty()) gflow->setAllBCs(choose_bc(h->params[0]->partA));
        else {
          for (auto bc : h->subHeads) {
            if (!h->heading.empty() || h->params.size()==1) { // A specific dimension must be choosen, the (one) parameter is the flag
              gflow->setBC(convert<int>(h->heading), choose_bc(h->params[0]->partA));
            }
            else throw BadStructure();
          }
        }
      }
      else good = false;
    }

    // Look for fill areas
    token = "Fill-area:";
    good = true;
    it = options.find(token);
    // Look for options
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // This is complicated enough that we give it it's own function. Fill area has its own options.
        fillArea(h);
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
    if (token=="HardSphere") return new HardSphere(gflow);
    else throw UnexpectedOption();
  }

  inline BCFlag FileParseCreator::choose_bc(string& token) const {
    if (token=="Wrap")      return BCFlag::WRAP;
    else if (token=="Repl") return BCFlag::REPL;
    else throw UnexpectedOption();
  }

  inline void FileParseCreator::fillArea(HeadNode *head) const {
    string token;
    bool good = true;
    std::multimap<string, HeadNode*>::iterator it;

    // --- Necessary data
    Bounds bounds; // We must have bounds
    int number = 0; // We must initialize a positive number of particles
    bool usePhi = false;
    RealType phi = 0;
    // Selection functions
    std::function<int(int)> select_type = nullptr;
    std::function<RealType(int)> select_sigma = nullptr;
    std::function<RealType(int, RealType)> select_mass = nullptr;
    std::function<void(RealType*, RealType*, RealType, RealType, int)> select_velocity = nullptr;
    // Data vector for lambda functions to reference
    vector<RealType> data;

    // --- Sort options
    std::map<string, HeadNode*> options;
    for (auto h : head->subHeads)
      options.insert(std::pair<string, HeadNode*>(h->heading, h));

    // --- Look for options

    token = "Bounds:";
    good = true;
    it = options.find(token);
    // Look for options
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        // The body of bounds contains the actual bounds
        if (h->subHeads.size()!=DIMENSIONS) throw BadStructure();
        // Set bounds
        for (int d=0; d<h->subHeads.size(); ++d) {
          if (h->subHeads[d]->params.size()!=2) throw BadStructure();
          bounds.min[d] = convert<float>( h->subHeads[d]->params[0]->partA );
          bounds.max[d] = convert<float>( h->subHeads[d]->params[1]->partA );
        }
      }
      else good = false;
    }

    token = "Number:";
    good = true;
    it = options.find(token);
    // Look for options
    for (; it!=options.end() && good; ++it) {
      if (it->first==token) {
        HeadNode *h = it->second;
        if (h->params.size()!=1) throw BadStructure();
        string a = h->params[0]->partA, b = h->params[0]->partB;
        if (!b.empty()) { // There is an option
          if (a=="Phi") {
            usePhi = true;
            phi = convert<RealType>(b);
          }
          else throw UnexpectedOption();
        }
        // Just read the number
        else number = convert<int>(a);
      }
      else good = false;
    }

    // FOR NOW 
    // All particles are of type 0
    select_type = [] (int) -> int {
      return 0;
    };
    // Constant radius - radius is 0.05
    select_sigma = [] (int) -> RealType {
      return 0.05;
    };
    // Constant density - density is 1
    select_mass = [] (int, RealType sigma) -> RealType {
      return 1.0 * sphere_volume(sigma);
    };
    // Select a velocity
    select_velocity = [&] (RealType *V, RealType *X, RealType sigma, RealType im, int type) -> void {
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
    if (select_type==nullptr || select_sigma==nullptr || select_mass==nullptr || select_velocity==nullptr) return;

    // --- We have found all the options. Fill the area.
    GFlow *filler = new GFlow;
    filler->setBounds(bounds);
    filler->setAllBCs(BCFlag::WRAP);
    filler->forceMaster = gflow->forceMaster; // Make sure the particles treat each other in the same way
    // Get the simdata
    SimData *simData = filler->simData;
    // --- Fill with particles
    RealType X[DIMENSIONS], V[DIMENSIONS];
    zeroVec(V);
    // If we are filling to a specified packing fraction
    if (usePhi) {
      RealType vol = 0, Vol = bounds.vol();
      int i=0; // A counter
      while (vol/Vol < phi) {
        // Select a position for the particle (random uniform)
        for (int d=0; d<DIMENSIONS; ++d)
          X[d] = drand48()*bounds.wd(d) + bounds.min[d];
        // Select other characteristics
        int type = select_type(i);
        RealType sigma = select_sigma(i);
        RealType im = 1./select_mass(i, sigma);
        // Add the particle
        simData->addParticle(X, V, sigma, im, type);
        // Increment volume and counter
        vol += sphere_volume(sigma);
        ++i;
      }
    }
    // If we are filling to a specified number
    else {
      for (int i=0; i<number; ++i) {
        // Select a position for the particle (random uniform)
        for (int d=0; d<DIMENSIONS; ++d)
          X[d] = drand48()*bounds.wd(d) + bounds.min[d];
        // Select other characteristics
        int type = select_type(i);
        RealType sigma = select_sigma(i);
        RealType im = 1./select_mass(i, sigma);
        // Add the particle
        simData->addParticle(X, V, sigma, im, type);
      }
    }
    // Initialize domain
    filler->domain->initialize();

    // --- Relax the simulation
    hs_relax(filler, 0.1); // 1) To make sure particles don't stop on top of one another
    relax(filler, 0.15);

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

    // Delete 
    delete filler;
  }

  inline void FileParseCreator::passComment(std::ifstream& fin, bool mline) {
    // Start right after "//" or "/*"
    char c;
    fin.get(c);
    while (!fin.eof()) {
      if (mline && c=='*') {
        fin.get(c);
        if (fin.eof()) return; // Check end of file
        if (c=='/') return; // End of the comment
      }
      else if (c=='\n' || c=='\r') {// Single line comments end with at newline
        fin.putback(c); // In case a newline signifies something for whoever called this function
        return;
      }
      // Get next character
      fin.get(c);
    }
  }

  inline bool FileParseCreator::passSpaces(std::ifstream& fin) {
    char c;
    // Get first character
    fin.get(c);
    // Loop
    while (!fin.eof()) {
      if (c==' ');
      else if (c=='\n' || c=='\r') return true;
      else { // Encountered a non-whitespace
        fin.putback(c);
        return false;
      }
      // Get next character
      fin.get(c);
    }
  }

  inline void FileParseCreator::passWhiteSpaces(std::ifstream& fin) {
    char c;
    // Get first character
    fin.get(c);
    // Loop
    while (!fin.eof()) {
      if (c==' ' || c=='\n' || c=='\r');
      else { // Encountered a non-whitespace
        fin.putback(c);
        return;
      }
      // Get next character
      fin.get(c);
    }
  }

  inline void FileParseCreator::getBody(std::ifstream& fin) {
    // Look for heads
    char c;
    bool end = false;
    while (!fin.eof() && !end) {
      passWhiteSpaces(fin);

      if (!fin.eof()) fin.get(c);
      else return;

      if (c=='}') // End of a body
        return;
      else if (c=='/') { // Could be the start of a comment
        checkComment(fin);
      }
      else {
        fin.putback(c);
        getHead(fin);
      }
    }
  }

  inline void FileParseCreator::getHead(std::ifstream& fin) {
    if (fin.eof()) return;
    // Create a new head node
    HeadNode *node = new HeadNode;

    // Get the heading
    string str;
    fin >> str;
    // Should end with a ":"
    node->heading = str;

    // Set node as the current head
    node->parent = currentHead;
    currentHead = node;

    // Look for parameters
    bool newLine = passSpaces(fin);
    bool body = false;

    // Get the parameters - adds them to the current head node
    if (!fin.eof() && !newLine)
      body = getParameters(fin);

    ++level;
    if (body && !fin.eof()) getBody(fin);
    --level;

    // Add this node to the parent node
    node->parent->subHeads.push_back(node);

    // Return to parent node
    currentHead = node->parent;
  }

  inline bool FileParseCreator::getParameters(std::ifstream& fin) {    
    bool more = true, body = false;
    while (more) getParam(fin, more, body);
    // Return whether we expect a body or not
    return body;
  }

  inline void FileParseCreator::getParam(std::ifstream& fin, bool& more, bool& body) {
    char c;
    bool end = false, a_part = true;
    fin.get(c);
    string a(""), b("");
    while (!fin.eof()) {
      if (c=='/') {
        checkComment(fin);
      }
      else if (c=='=') {
        a_part = false;
      }
      else if (c==',') {
        more = true;
        break;
      }
      else if (c=='{') {
        more = false;
        body = true;
        break;
      }
      else if (c=='\n' || c=='\r') {
        body = false;
        more = false;
        break;
      }
      else if (c==' ');
      else if (c=='\n' || c=='\r') {
        more = false;
        body = false;
        break;
      }
      else {
        if (a_part) a.push_back(c);
        else        b.push_back(c);
      }

      // Get next character
      if (!fin.eof()) fin.get(c);
    }
    // Set param
    if (!a.empty() || !b.empty()) currentHead->params.push_back(new ParamNode(a, b));
  }

  inline void FileParseCreator::checkComment(std::ifstream& fin) {
    char c;
    fin.get(c);
    // Make sure we are not at eof
    if (fin.eof()) return;
    // Check what kind of comment this is
    if (c=='/')      passComment(fin, false);
    else if (c=='*') passComment(fin, true);
    else             throw UnexpectedToken();
  }

}