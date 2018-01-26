#include "FileParser.hpp"
#include "Creator.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"

namespace GFlow {

  FileParser::FileParser(int argc, char** argv) : creator(new Creator), status(nullptr), randomSeed(0), argc(argc), argv(argv), parser(nullptr) {};

  FileParser::~FileParser() {
    if (creator) delete creator;
  }
  
  void FileParser::parse(string filename, SimData *& simData, Integrator *& integrator, unsigned seed) {
    // Timing
    if (dataRecord) dataRecord->startTiming();

    // Open file
    ifstream fin(filename);
    if (fin.fail()) {
      if (status) status->writeError("File ["+filename+"] cannot be opened by file parser.");
      throw FileDoesNotExist(filename);
    }
    if (status) status->writeMessage("File ["+filename+"] opened by file parser. Will attempt to construct configuration from file.");
    
    // Seed random number generator
    if (seed==0) seed_rand();
    else seed_rand(seed);

    // Data to save
    Bounds bounds(NullBounds);
    vector<Particle> particles;
    vector<Region> regions;
    vector<Wall> walls;
    vec2 gravity;
    bool wrapX(false), wrapY(false);

    // Integrator data to save
    RealType dt(-1), minDt(-1), maxDt(-1);
    bool adjustDt(true), adjustDelay(true);

    // For user defined variables and options
    std::map<string, string> variables;
    std::map<string, string> options;
    // Define true as 1 and false as 0
    variables.emplace(Bool_True,  "1");
    variables.emplace(Bool_False, "0");

    // First look for options
    find_options(filename, options);

    // Check if the options have been passed in through the command line
    bool createdParser = false;
    if (parser==nullptr) {
      parser = new ArgParse;
      createdParser = true;
      try {
	parser->set(argc, argv);
      }
      catch (ArgParse::IllegalToken token) {
	std::cerr << "Illegal Argument: " << token.c << ". Exiting.\n";
	exit(1);
      }
    }
    for (const auto &p : options) {
      string value;
      parser->get(p.first, value);
      // If a value is supplied from the command line
      if (value!="") variables.emplace(p.first, value);
      // Else, use the default value
      else variables.emplace(p.first, p.second);
    }
    if (createdParser) delete parser;

    // Helper data
    const int max_comment_size = 512;
    char comment[max_comment_size];
    string tok, opt;
    // Get the first token
    fin >> tok;
    while (!fin.eof()) {      
      // (// ... \ncomment) // 
      if (tok.at(0)==Comment_Tok && tok.at(1)==Comment_Tok) fin.getline(comment, max_comment_size);
      else if (tok.at(0)=='/' && tok.at(1)=='*') {
	char f, s;
	fin.get(f); fin.get(s);
	while (!fin.eof()) {
	  if (f=='*' && s=='/') break;
	  f = s;
	  fin.get(s);
	}
      }
      else if (tok==Variable_Tok) { // let [name] [value]
	string name, val;
	fin >> name;
	getValue(fin, val, variables);
	variables.emplace(name, val);
      }
      else if (tok==Option_Tok) {
	string name, val;
	fin >> name >> val;
	// We already got all the options, so we can ignore it here
      }
      else if (tok==Dt_Tok) {
	getValue(fin, dt, variables);
      }
      else if (tok==MinDt_Tok) {
	getValue(fin, minDt, variables);
      }
      else if (tok==MaxDt_Tok) {
	getValue(fin, maxDt, variables);
      }
      else if (tok==AdjustDt_Tok) {
	getValue(fin, adjustDt, variables);
      }
      else if (tok==AdjustDelay_Tok) {
	getValue(fin, adjustDelay, variables);
      }
      else if (tok==WrapX_Tok) { // wrapX [true/false]
	bool wrap;
	getValue(fin, wrap, variables);
	wrapX = wrap;
      }
      else if (tok==WrapY_Tok) { // wrapY [true/false]
	bool wrap;
	getValue(fin, wrap, variables);
        wrapY = wrap;
      }
      else if (tok==Bounds_Tok) { // bounds [lf] [rt] [bt] [tp]
	RealType lf, rg, bt, tp;
	getValue(fin, lf, variables);
	getValue(fin, rg, variables);
	getValue(fin, bt, variables);
	getValue(fin, tp, variables);
	bounds = Bounds(lf, rg, bt, tp);
      }
      else if (tok==Wall_Tok) make_wall(fin, walls, variables); // wall
      else if (tok==Particle_Tok) make_particle(fin, particles, variables); // particle
      else if (tok==Gravity_Tok) { // gravity [ax] [ay]
	RealType x,y;
	getValue(fin, x, variables);
	getValue(fin, y, variables);
	gravity = vec2(x,y);
      }
      
      // region 
      // [left] [right] [bottom] [top]
      // (options) N [number] | disp [dispersion] | disp_type [dispersion type] | it [interaction] | ...
      // end
      else if (tok==Region_Tok) make_region(fin, regions, variables);
      // Unrecognized
      else throw UnrecognizedToken(tok);      

      // Get the next character
      fin >> tok;
    }
    
    // Create simulation data
    if (bounds==NullBounds) {
      simData = nullptr;
      integrator = nullptr;
      return;
    }
    
    // Set simulation bounds
    simData = new SimData(bounds, bounds);
    // Set "gravity"
    if (gravity!=Zero) simData->addExternalForce(new ConstantAcceleration(gravity));
    // Set wrapping
    simData->setWrapX(wrapX);
    simData->setWrapY(wrapY);
    // Add walls
    simData->addWall(walls);
    // Add individual particles
    simData->reserve(particles.size());
    for (auto& p : particles) simData->addParticle(p);
    // Add particles from regions
    for (auto& r : regions) {
      // Set bounds to the entire region if the bounds are unspecified
      if (r.bounds==NullBounds) r.bounds = bounds;
      creator->createRegion(r, simData);
    }

    // Create integrator
    VelocityVerletIntegrator* integ = new VelocityVerletIntegrator(simData);
    integrator = integ;
    if (dt>0)    integ->setDt(dt);
    integ->setAdjustTimeStep(adjustDt);
    integ->setAdjustUpdateDelay(adjustDelay);
    if (minDt>0) integ->setMinTimeStep(minDt);
    if (maxDt>0) integ->setMaxTimeStep(maxDt);

    // Make sure there is no serious overlap
    StandardSectorization remover(simData);
    remover.removeOverlapping();

    // Timing
    if (dataRecord) {
      dataRecord->endTiming();
      RealType time = dataRecord->getElapsedTime();
      dataRecord->setSetupTime(time);
    }
  }

  SimData* FileParser::loadLegacyFromFile(string filename) {
    // Timing
    if (dataRecord) dataRecord->startTiming();

    // Open stream, check if failed
    ifstream fin(filename);
    if (fin.fail()) {
      if (status) status->writeError("File ["+filename+"] cannot be opened by file parser.");
      throw FileDoesNotExist(filename);
    }
    if (status) status->writeMessage("File ["+filename+"] opened by file parser. Will attempt to load saved state from file.");

    // Simulation data to look for
    RealType left, right, bottom, top;
    vec2 gravity;
    vector<RealType> radii;
    vector<vec2> wallLeft, wallRight, positions;
    vector<int> interactions;

    // Get simulation bounds
    fin >> left >> right >> bottom >> top;

    Bounds bounds(left, right, bottom, top);
    SimData *simData = new SimData(bounds, bounds);

    // Get gravity
    fin >> gravity;
    if (gravity!=Zero) 
      simData->addExternalForce(new ConstantAcceleration(gravity));
    
    // Get walls
    fin >> wallLeft;
    fin >> wallRight;

    // Get radii and positions
    fin >> radii;
    fin >> positions;
    fin >> interactions;
    fin.close();

    // Create walls and particles
    int size = positions.size(), rsize = radii.size(), wsize = min(wallLeft.size(), wallRight.size()), isize = interactions.size();
    for (int i=0; i<wsize; ++i)
      simData->addWall(Wall(wallLeft[i], wallRight[i]));
    // Reserve space
    simData->reserveAdditional(size);
    for (int i=0; i<size; ++i) {
      Particle p(positions.at(i), radii.at(i%rsize));
      p.interaction = interactions.at(i%isize);
      simData->addParticle(p);
    }

    // Timing
    if (dataRecord) {
      dataRecord->endTiming();
      RealType time = dataRecord->getElapsedTime();
      dataRecord->setSetupTime(time);
    }

    // Return the sim data object
    return simData;
  }

  void FileParser::loadFromFile(string filename, SimData *& simData, Integrator *& integrator) {
    // Timing
    if (dataRecord) dataRecord->startTiming();
    // Open the file
    ifstream fin(filename);
    if (fin.fail()) {
      if (status) status->writeError("File ["+filename+"] cannot be opened by file parser when attempting to load arrangement.");
      throw FileDoesNotExist(filename);
    }
    
    // Data to save
    Bounds bounds(NullBounds);
    vector<Particle> particles;
    vector<Wall> walls;
    vector<ExternalForce*> externalForces;
    vector<int> fixedParticles;
    vector<pair<vec2, vector<int> > > drivenParticles;
    vector<pair<vec2, vector<int> > > cvParticles;
    bool wrapX(false), wrapY(false);

    // For getting comments
    const int max_comment_size = 512;
    char comment[max_comment_size];

    // Get tokens and data
    string tok;
    fin >> tok;
    while (!fin.eof()) {
      // Comments
      if (tok.at(0)=='#') { // Comment
	fin.getline(comment, max_comment_size);
      }
      // Bounds
      if (tok=="B") {
	RealType l, r, b, t;
	fin >> l >> r >> b >> t;
	bounds = Bounds(l,r,b,t);
      }
      // Whether we should wrap in the x direction (0/1)
      else if (tok=="wx") {
	fin >> wrapX;
      }
      // Whether we should wrap in the y direction (0/1)
      else if (tok=="wy") {
	fin >> wrapY;
      }
      // Add an external force
      else if (tok=="ef") {
	fin >> tok; // Get the type of external force
	if (tok=="ca") {
	  vec2 acc;
	  fin >> acc;
	  externalForces.push_back(new ConstantAcceleration(acc));
	}
	if (tok=="vd") {
	  RealType vis;
	  fin >> vis;
	  externalForces.push_back(new ViscousDrag(vis));
	}
      }
      // Create a wall
      else if (tok=="W") {
	vec2 left, right;
	RealType rp, ds, cf;
	// Get the wall data
	fin >> left >> right >> rp >> ds >> cf;
	Wall w(left, right);
	w.repulsion   = rp;
	w.dissipation = ds;
	w.coeff       = cf;
	// Push back the wall
	walls.push_back(w);
      }
      // Create a particle
      else if (tok=="P") {
	vec2 pos, velocity;
	RealType th, om, sg, rp, ds, cf;
	int it;
	fin >> pos >> velocity >> th >> om >> sg >> rp >> ds >> cf >> it;
	// Set the particle data
	Particle P(pos, sg);
	P.velocity    = velocity;
	P.theta       = th;
	P.omega       = om;
	P.repulsion   = rp;
	P.dissipation = ds;
	P.coeff       = cf;
	P.interaction = it;
	// Push back the particle
	particles.push_back(P);
      }
      // Fixed particles
      else if (tok=="FP") {
	vector<int> particleIDs;
	fin >> particleIDs;
	fixedParticles.insert(fixedParticles.end(), particleIDs.begin(), particleIDs.end());
      }
      // External force applied to particles
      else if (tok=="EFP") {
	vector<int> particleIDs;
	vec2 df;
	// Read in the ids of the particles and the df
	fin >> particleIDs;
	fin >> df;
	drivenParticles.push_back(pair<vec2, vector<int> >(df, particleIDs));
      }
      // Constant velocity particles
      else if (tok=="CV") {
	vector<int> particleIDs;
	vec2 vel;
	fin >> particleIDs;
	fin >> vel;
	cvParticles.push_back(pair<vec2, vector<int> >(vel, particleIDs));
      }
      // Get the next token
      fin >> tok;
    }
    fin.close();

    // Create the simulation data
    if (bounds==NullBounds) {
      simData = nullptr;
      integrator = nullptr;
      return;
    }
    // Set simulation bounds
    simData = new SimData(bounds, bounds);
    // Set forces
    for (auto ef : externalForces)
      simData->addExternalForce(ef);
    // Set wrapping
    simData->setWrapX(wrapX);
    simData->setWrapY(wrapY);
    // Add walls
    for (auto& w : walls) simData->addWall(w);
    // Add individual particles
    simData->reserve(particles.size());
    simData->addParticle(particles);
    /** Add characteristics as necessary **/
    // Add fixed particles
    for (auto id : fixedParticles)
      simData->addCharacteristic(id, new Fixed(vec2(simData->getPxPtr()[id], simData->getPyPtr()[id])));
    // Add driven particles
    for (auto pr : drivenParticles) {
      auto df = pr.first;
      for (auto id : pr.second)
	simData->addCharacteristic(id, new ApplyForce(Zero, df));
    }
    // Add constant velocity particles
    for (auto pr : cvParticles) {
      auto v = pr.first;
      for (auto id : pr.second) 
	simData->addCharacteristic(id, new ConstantVelocity(vec2(simData->getPxPtr() [id], simData->getPyPtr() [id]), v));
    }
    
    // Create integrator
    integrator = new VelocityVerletIntegrator(simData);

    // Timing
    if (dataRecord) {
      dataRecord->endTiming();
      RealType time = dataRecord->getElapsedTime();
      dataRecord->setSetupTime(time);
    }
  }

  void FileParser::saveToFile(SimData* simData, string filename) {
    ofstream fout(filename);
    if (fout.fail()) {
      if (status) status->writeError("File ["+filename+"] cannot be opened by file parser when attempting to save arrangement.");
      throw FileDoesNotExist(filename);
    }
    // Write number of balls and walls in a comment
    fout << "# There are " << simData->domain_size << " particles and " << simData->walls.size() << ".\n";
    // Write bounds
    Bounds simBounds = simData->simBounds;
    fout << "B " << simBounds.left << " " << simBounds.right << " " << simBounds.bottom << " " << simBounds.top << endl;
    // Write wrapping
    fout << "wx " << (simData->wrapX ? "1" : "0") << "\nwy " << (simData->wrapY ? "1" : "0") << "\n";
    // Write walls
    for (const auto& w : simData->walls) {
      fout << "W " << w.left << " " << w.getRight() << " " << w.repulsion << " " << w.dissipation << " " << w.coeff << "\n";
    }
    // Write external forces
    for (const auto& f : simData->externalForces) {
      auto ca = reinterpret_cast<ConstantAcceleration*>(f);
      if (ca) fout << "ef ca " << ca->getAcceleration() << "\n";
      // Allow for other forces here
    }
    // Write particles
    for (const auto& p : simData->getParticles()) {
      fout << "P " << p.position << " " << p.velocity << " " << p.theta << " " << p.omega << " " << p.sigma << " " << p.repulsion << " " << p.dissipation << " " << p.coeff << " " << p.interaction << "\n";
    }
    // Nothing else to write. Close stream.
    fout.close();
  }
  
  // Make a region data structure
  inline void FileParser::make_region(std::ifstream& fin, vector<Region>& regions, const std::map<string,string>& variables) {    
    // region
    // bounds [left] [right] [bottom] [top] | number [number] | disp [dispersion] | disp_type [dispersion type] | it [interaction] | ...
    // ...
    // end

    // For recording data
    string tok("");
    RealType sigma(-1), dispersion(0), phi(-1), coeff(-1), repulsion(-1), dissipation(-1), velocity(-1), vy(-1);
    Bounds bounds = NullBounds;
    int number(-1), interaction(-1);
    Region region;
    
    // Choice of how we choose velocities
    bool randomVChoice[3]; 
    for(int i=0; i<3; ++i) randomVChoice[i] = false;

    // For getting comments
    const int max_comment_size = 512;
    char comment[max_comment_size];

    // Loop
    bool end = false;
    fin >> tok;
    while (!fin.eof() && !end) {
      if (tok.at(0)==Comment_Tok) {
        fin.getline(comment, max_comment_size);
      }
      else if (tok==End_Tok) {
	end = true;
	continue;
      }
      else if (tok==Number_Tok) getValue(fin, number, variables);
      else if (tok==Phi_Tok) getValue(fin, phi, variables);
      else if (tok==Sigma_Tok) getValue(fin, sigma, variables);
      else if (tok==Dispersion_Tok) getValue(fin, dispersion, variables);
      else if (tok==Normal_Velocity_Tok) {
	getValue(fin, velocity, variables);
	randomVChoice[0] = true;
      }
      else if (tok==Normal_KE_Tok) {
	getValue(fin, velocity, variables);
	randomVChoice[1] = true;
      }
      else if (tok==Velocity_Tok) {
	getValue(fin, velocity, variables);
	getValue(fin, vy, variables);
	randomVChoice[2] = true;
      }
      else if (tok==Repulsion_Tok) getValue(fin, repulsion, variables);
      else if (tok==Dissipation_Tok) getValue(fin, dissipation, variables);
      else if (tok==Coeff_Tok) getValue(fin, coeff, variables);
      else if (tok==Interaction_Tok) getValue(fin, interaction, variables);
      else if (tok==Bounds_Tok) {
	RealType left, right, bottom, top;
	getValue(fin, left, variables);
	getValue(fin, right, variables);
	getValue(fin, bottom, variables);
	getValue(fin, top, variables);
	bounds = Bounds(left, right, bottom, top);
      }
      else throw UnrecognizedToken(tok);
      
      // Get next token
      fin >> tok;
    }

    // Sigma function
    if (number!=-1) {
      region.sigma = new Fixed_Number_Uniform_Radii(number, sigma, dispersion);
    }
    else if (phi!=-1) {
      region.sigma = new Fixed_Phi_Uniform_Radii(phi, sigma, dispersion);
    }
    // Dissipation
    if (0<=dissipation) {
      region.dissipation = new Uniform_Random_Dissipation(dissipation);
    }
    // Velocity
    if (0<=velocity) {
      if (randomVChoice[0])      region.velocity = new Normal_Random_Velocity(velocity);
      else if (randomVChoice[1]) region.velocity = new Normal_Random_KE(velocity);
      else if (randomVChoice[2]) region.velocity = new Constant_Velocity(velocity, vy);
    }
    // Repulsion
    if (0<=repulsion) {
      region.repulsion = new Uniform_Random_Repulsion(repulsion);
    }
    // Coefficient of friction
    if (0<=coeff) {
      region.coeff = new Uniform_Random_Coeff(coeff);
    }
    // Interaction
    if (-1<interaction) {
      region.interaction = new Homogeneous_Interaction(interaction);
    }
    
    // Bounds will be NullBounds (whole space) if no bounds were specified
    region.bounds = bounds;

    // Add the region
    regions.push_back(region);
  }

  void FileParser::make_wall(std::ifstream& fin, vector<Wall>& walls, const std::map<string,string>& variables) {

    // wall
    // pos [lx] [ly] [rx] [ry] | repulsion [repulsion] | dissipation [dissipation] | coeff [coeff]
    // ...
    // end    

    // For recording data
    string tok("");
    RealType lx(0), ly(0), rx(0), ry(0);
    RealType repulsion(default_wall_repulsion), dissipation(default_wall_dissipation), coeff(default_wall_coeff);

    // For getting comments
    const int max_comment_size = 512;
    char comment[max_comment_size];

    // Loop
    bool end = false;
    fin >> tok;
    while (!fin.eof() && !end) {
      if (tok.at(0)==Comment_Tok) {
        fin.getline(comment, max_comment_size);
      }
      else if (tok==End_Tok) {
	end = true;
	continue;
      }
      else if (tok==Position_Tok) {
	getValue(fin, lx, variables);
	getValue(fin, ly, variables);
	getValue(fin, rx, variables);
	getValue(fin, ry, variables);
      }
      else if (tok==Repulsion_Tok) getValue(fin, repulsion, variables);
      else if (tok==Dissipation_Tok) getValue(fin, dissipation, variables);
      else if (tok==Coeff_Tok) getValue(fin, coeff, variables);
      else throw UnrecognizedToken(tok);
      // Get next token
      fin >> tok;
    }

    // Make the wall from the gathered values
    if (!(lx==0 && ly==0 && rx==0 && ry==0)) {
      Wall w(lx, ly, rx, ry);
      w.repulsion = repulsion;
      w.dissipation = dissipation;
      w.coeff = coeff;
      walls.push_back(w);
    }
  }

  void FileParser::make_particle(std::ifstream& fin, vector<Particle>& particles, const std::map<string, string>& variables) {
    // wall
    // pos [lx] [ly] [rx] [ry] | repulsion [repulsion] | dissipation [dissipation] | coeff [coeff]
    // ...
    // end
    
    // For recording data
    string tok("");
    RealType X(0), Y(0), R(-1), vx(0), vy(0), omega(0);
    bool gotPos = false;
    RealType repulsion(default_particle_repulsion), dissipation(default_particle_dissipation), coeff(default_particle_coeff);
    
    // For getting comments
    const int max_comment_size = 512;
    char comment[max_comment_size];
    
    // Loop
    bool end = false;
    fin >> tok;
    while (!fin.eof() && !end) {
      if (tok.at(0)==Comment_Tok) {
	fin.getline(comment, max_comment_size);
      }
      else if (tok==End_Tok) {
	end = true;
	continue;
      }
      else if (tok==Position_Tok) {
	getValue(fin, X, variables);
	getValue(fin, Y, variables);
	gotPos = true;
      }
      else if (tok==Sigma_Tok) getValue(fin, R, variables);
      else if (tok==Repulsion_Tok) getValue(fin, repulsion, variables);
      else if (tok==Dissipation_Tok) getValue(fin, dissipation, variables);
      else if (tok==Coeff_Tok) getValue(fin, coeff, variables);
      else if (tok==Velocity_Tok) {
	getValue(fin, vx, variables);
	getValue(fin, vy, variables);
      }
      else if (tok==Omega_Tok) getValue(fin, omega, variables);
      else throw UnrecognizedToken(tok);
      
      fin >> tok;
    }
    
    // Create and add the particle to the list
    if (gotPos && 0<R) {
      Particle P(X,Y,R);
      P.velocity = vec2(vx, vy);
      P.omega = omega;
      P.repulsion = repulsion;
      P.dissipation = dissipation;
      P.coeff = coeff;
      particles.push_back(P);
    }
  }

  inline void FileParser::find_options(string filename, std::map<string, string>& options) {
    std::ifstream fin(filename);
    if (fin.fail()) {
      std::cerr << "Failure.\n";
      exit(0);
    }
    // Get all options
    string tok;
    fin >> tok;
    while (!fin.eof()) {
      if (tok==Option_Tok) {
	string value;
	fin >> tok >> value;
	options.emplace(tok, value);
      }
      fin >> tok;
    }
  }
  
  inline void FileParser::seed_rand() {
    randomSeed = std::chrono::system_clock::now().time_since_epoch().count();
    srand48( randomSeed );
    srand  ( randomSeed );
    generator = std::mt19937(randomSeed);
  }
  
  inline void FileParser::seed_rand(unsigned seed) {
    randomSeed = seed;
    srand48( randomSeed );
    srand  ( randomSeed );
    generator =std::mt19937(randomSeed);
  }

}
