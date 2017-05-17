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
    fin >> c;

    while (!fin.eof()) {      
      // (comment) // 
      if (c=='/') {
	fin.getline(comment, max_comment_size);
	fin.get(c);
	continue;
      }
      
      // Otherwise
      if (isalpha(c)) {
	fin.putback(c);
	fin >> tok;
      }
      else { // Error
	
      }
      
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
	else { // Error
	  
	}
      }
      
      // wall [lx] [ly] [rx] [ry]
      else if (tok=="wall") {
	RealType lx, ly, rx, ry;
	fin >> lx >> ly >> rx >> ry;
	walls.push_back(Wall(lx, ly, rx, ry));
      }
      
      // region 
      // [left] [right] [bottom] [top]
      // (options) N [number] | disp [dispersion] | disp_type [dispersion type] | it [interaction] | ...
      // end
      else if (tok=="region") make_region(fin);
      
      // Get the next character
      fin >> c;
    }

    // Create simulation data
    if (simBounds==NullBounds) return nullptr;
    simData->setBounds(simBounds);
    simData->addExternalForce(new ConstantAcceleration(gravity));
    simData->setWrapX(wrapX);
    simData->setWrapY(wrapY);
    for (auto& w : walls) {
      simData->addWall(w);
    }
    simData->reserve(particles.size());
    for (auto& p : particles) {
      simData->addParticle(p);
    }
    for (auto& r : regions) {
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
    RealType r_value(0), sigma(-1), dispersion(0), phi(-1);
    int i_value(0), number(-1);
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
      else if (tok=="sigma" || tok=="radius") {
	fin >> sigma;
      }
      else if (tok=="dispersion") {
	fin >> dispersion;
      }
      else if (tok=="phi") {
	fin >> phi;
      }

      // Get next token
      fin >> tok;
    }

    if (number!=-1) {
      if (region.sigma) delete region.sigma;
      region.sigma = new Fixed_Number_Uniform_Radii(number, sigma, dispersion);
    }
    else if (phi!=-1) {
      if (region.sigma) delete region.sigma;
      region.sigma = new Fixed_Phi_Uniform_Radii(phi, sigma, dispersion);
    }

    // Feed the Region to the creator
    regions.push_back(region);
  }

}
