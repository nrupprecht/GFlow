#include "FileParser.hpp"
#include "Creator.hpp"

namespace GFlow {

  FileParser::FileParser() : simData(nullptr), creator(new Creator), randomSeed(0) {
    // Default values
    gravity   = Zero;
    wrapX     = false;
    wrapY     = false;
    simBounds = NullBounds;
  };

  FileParser::~FileParser() {
    if (creator) delete creator;
  }
  
  SimData* FileParser::parse(string filename, unsigned seed) {
    ifstream fin(filename);
    if (fin.fail()) {
      throw FileDoesNotExist(filename);
    }

    // Seed random number generator
    if (seed==0) seed_rand();
    else seed_rand(seed);

    // Helper data
    const int max_comment_size = 512;
    char comment[max_comment_size];
    string tok, opt;
    // Get the first token
    fin >> tok;
    while (!fin.eof()) {      
      // (comment) // 
      if (tok.at(0)==Comment_Tok) {
	fin.getline(comment, max_comment_size);
	fin >> tok;
	continue;
      }
      
      // Otherwise
      if (!isalpha(tok.at(0))) throw UnrecognizedToken(tok);

      // wrapX [true/false]
      if (tok==WrapX_Tok || tok==WrapY_Tok) {
	fin >> opt;
	if (opt==Bool_True) {
	  if (tok==WrapX_Tok) wrapX = true;
	  else              wrapY = true;
	}
	else if (opt==Bool_False) {
	  if (tok==WrapX_Tok) wrapX = false;
          else              wrapY = false;
	}
	else return nullptr;
      }

      // bounds [lf] [rt] [bt] [tp]
      else if (tok==Bounds_Tok) {
	RealType lf, rg, bt, tp;
	fin >> lf >> rg >> bt >> tp;
	simBounds = Bounds(lf, rg, bt, tp);
      }
      
      // wall
      else if (tok==Wall_Tok) make_wall(fin);

      // particle
      else if (tok==Particle_Tok) make_particle(fin);

      // gravity [ax] [ay]
      else if (tok==Gravity_Tok) {
	RealType x,y;
	fin >> x >> y;
	gravity = vec2(x,y);
      }
      
      // region 
      // [left] [right] [bottom] [top]
      // (options) N [number] | disp [dispersion] | disp_type [dispersion type] | it [interaction] | ...
      // end
      else if (tok==Region_Tok) make_region(fin);
      
      // Get the next character
      fin >> tok;
    }

    // Create simulation data
    if (simBounds==NullBounds) return nullptr;
    // Set simulation bounds
    simData = new SimData(simBounds, simBounds);
    // Set "gravity"
    if (gravity!=Zero) simData->addExternalForce(new ConstantAcceleration(gravity));
    // Set wrapping
    simData->setWrapX(wrapX);
    simData->setWrapY(wrapY);
    // Add walls
    for (auto& w : walls) simData->addWall(w);
    // Add individual particles
    simData->reserve(particles.size());
    for (auto& p : particles) simData->addParticle(p);
    // Add particles from regions
    for (auto& r : regions) {
      // Set bounds to the entire region if the bounds are unspecified
      if (r.bounds==NullBounds) r.bounds = simBounds;
      creator->createRegion(r, simData);
    }

    return simData;
  }

  SimData* FileParser::loadFromFile(string filename) {
    // Open stream, check if failed
    ifstream fin(filename);
    if (fin.fail()) {
      throw FileDoesNotExist(filename);
    }

    // Data to look for
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
    for (int i=0; i<size; ++i) {
      Particle p(positions.at(i), radii.at(i%rsize));
      p.interaction = interactions.at(i%isize);
      simData->addParticle(p);
    }
    // Return the sim data object
    return simData;
  }

  void FileParser::saveToFile(SimData* simData, string filename) {
    
  }
  
  // Make a region data structure and store it
  inline void FileParser::make_region(std::ifstream& fin) {    
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
      else if (tok==Number_Tok) {
	fin >> number;
      }
      else if (tok==Phi_Tok) {
	fin >> phi;
      }
      else if (tok==Sigma_Tok) {
	fin >> sigma;
      }
      else if (tok==Dispersion_Tok) {
	fin >> dispersion;
      }
      else if (tok==Normal_Velocity_Tok) {
	fin >> velocity;
	randomVChoice[0] = true;
      }
      else if (tok==Normal_KE_Tok) {
	fin >> velocity;
	randomVChoice[1] = true;
      }
      else if (tok==Velocity_Tok) {
	fin >> velocity >> vy;
	randomVChoice[2] = true;
      }
      else if (tok==Repulsion_Tok) {
	fin >> repulsion;
      }
      else if (tok==Dissipation_Tok) {
	fin >> dissipation;
      }
      else if (tok==Coeff_Tok) {
	fin >> coeff;
      }
      else if (tok==Interaction_Tok) {
	fin >> interaction;
      }
      else if (tok==Bounds_Tok) {
	RealType left, right, bottom, top;
	fin >> left >> right >> bottom >> top;
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

    // Feed the Region to the creator
    regions.push_back(region);
  }

  inline void FileParser::make_wall(std::ifstream& fin) {

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
	fin >> lx >> ly >> rx >> ry;
      }
      else if (tok==Repulsion_Tok) {
	fin >> repulsion;
      }
      else if (tok==Dissipation_Tok) {
	fin >> dissipation;
      }
      else if (tok==Coeff_Tok) {
	fin >> coeff;
      }
      else throw UnrecognizedToken(tok);

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

  void FileParser::make_particle(std::ifstream& fin) {
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
	fin >> X >> Y;
	gotPos = true;
      }
      else if (tok==Sigma_Tok) {
	fin >> R;
      }
      else if (tok==Repulsion_Tok) {
	fin >> repulsion;
      }
      else if (tok==Dissipation_Tok) {
	fin >> dissipation;
      }
      else if (tok==Coeff_Tok) {
	fin >> coeff;
      }
      else if (tok==Velocity_Tok) {
	fin >> vx >> vy;
      }
      else if (tok==Omega_Tok) {
	fin >> omega;
      }
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
