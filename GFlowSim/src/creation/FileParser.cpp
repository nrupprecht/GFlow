#include "FileParser.hpp"
#include "Creator.hpp"

namespace GFlow {

  FileParser::FileParser() : creator(new Creator), simData(nullptr) {
    // Default values
    gravity   = Zero;
    wrapX     = false;
    wrapY     = false;
    simBounds = NullBounds;
  };

  FileParser::~FileParser() {
    if (creator) delete creator;
  }
  
  SimData* FileParser::parse(string filename) {
    ifstream fin(filename);
    if (fin.fail()) {
      return nullptr;
    }

    // Helper data
    const int max_comment_size = 512;
    char c;
    char comment[max_comment_size];
    string tok, opt;
    // Get the first token
    fin >> tok;
    while (!fin.eof()) {      
      // (comment) // 
      if (tok.at(0)=='/') {
	fin.getline(comment, max_comment_size);
	fin >> tok;
	continue;
      }
      
      // Otherwise
      if (!isalpha(tok.at(0))) return nullptr; // Error

      // wrapX [true/false]
      if (tok=="wrapX" || tok=="wrapY") {
	fin >> opt;
	if (opt=="true") {
	  if (tok=="wrapX") wrapX = true;
	  else              wrapY = true;
	}
	else if (opt=="false") {
	  if (tok=="wrapX") wrapX = false;
          else              wrapY = false;
	}
	else return nullptr;
      }

      // bounds [lf] [rt] [bt] [tp]
      else if (tok=="bounds") {
	RealType lf, rg, bt, tp;
	fin >> lf >> rg >> bt >> tp;
	simBounds = Bounds(lf, rg, bt, tp);
      }
      
      // wall [lx] [ly] [rx] [ry]
      else if (tok=="wall") {
	RealType lx, ly, rx, ry;
	fin >> lx >> ly >> rx >> ry;
	walls.push_back(Wall(lx, ly, rx, ry));
      }

      else if (tok=="gravity") {
	RealType x,y;
	fin >> x >> y;
	gravity = vec2(x,y);
      }
      
      // region 
      // [left] [right] [bottom] [top]
      // (options) N [number] | disp [dispersion] | disp_type [dispersion type] | it [interaction] | ...
      // end
      else if (tok=="region") make_region(fin);
      
      // Get the next character
      fin >> tok;
    }

    // Create simulation data
    if (simBounds==NullBounds) return nullptr;
    // Set simulation bounds
    simData = new SimData(simBounds, simBounds);
    // Set "gravity"
    if (gravity!=Zero) 
      simData->addExternalForce(new ConstantAcceleration(gravity));
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
  
  inline void FileParser::make_region(std::ifstream& fin) {
    // Make a region data structure and store it
    
    // region
    // [left] [right] [bottom] [top]
    // (options) number [number] | disp [dispersion] | disp_type [dispersion type] | it [interaction] | ...
    // end

    string tok("");
    RealType sigma(-1), dispersion(0), phi(-1), coeff(-1), repulsion(-1), dissipation(-1), velocity(-1);
    Bounds bounds = NullBounds;
    int number(-1), interaction(-1);
    Region region;

    bool end = false;
    fin >> tok;
    while (!fin.eof() && !end) {
      if (tok=="end") {
	end = true;
	continue;
      }
      else if (tok=="number") {
	fin >> number;
      }
      else if (tok=="phi") {
	fin >> phi;
      }
      else if (tok=="sigma" || tok=="radius") {
	fin >> sigma;
      }
      else if (tok=="dispersion") {
	fin >> dispersion;
      }
      else if (tok=="velocity") {
	fin >> velocity;
      }
      else if (tok=="repulsion") {
	fin >> repulsion;
      }
      else if (tok=="dissipation") {
	fin >> dissipation;
      }
      else if (tok=="coeff") {
	fin >> coeff;
      }
      else if (tok=="interaction") {
	fin >> interaction;
      }
      else if (tok=="bounds") {
	RealType left, right, bottom, top;
	fin >> left >> right >> bottom >> top;
	bounds = Bounds(left, right, bottom, top);
      }
      else throw false; // Invalid option

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
      region.velocity = new Normal_Random_Velocity(velocity);
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

}
