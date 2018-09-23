/*
 * Author: Nathaniel Rupprecht
 * Start Data: June 15, 2018
 *
 */

#include "../../include/ArgParse.hpp"
#include "../control/Sectorization.hpp"

using namespace GFlow;

int main(int argc, char** argv) {
  // Variables
  int trials = 1000;
  RealType minPhi = 0;
  RealType maxPhi = 0.1;
  RealType numPhi = 10;
  RealType sigma = 0.01;
  RealType size = 1.;
  
  // Read some simulation options
  ArgParse parser(argc, argv);
  parser.get("trials", trials);
  parser.get("minPhi", minPhi);
  parser.get("maxPhi", maxPhi);
  parser.get("numPhi", numPhi);
  parser.get("size", size);

  // Default bounds
  Bounds bounds(0., size, 0., size);

  // Area of one particle
  RealType area = PI*sqr(sigma);
  // D(Phi)
  RealType dPhi = (maxPhi-minPhi)/numPhi;
  // Store data
  vector<pair<RealType, RealType> > data;
  // Go through all density levels
  for (RealType phiTr = 0; phiTr<numPhi; ++phiTr) {
    RealType phi = phiTr*dPhi + minPhi;
    int number = sqr(size)*phi/area;

    int success = 0;
    for (int tr = 0; tr < trials; ++tr) {
      // Set up simulation data
      SimData simData(bounds, bounds);
      simData.setWrapX(true);
      simData.setWrapY(true);
      
      // Add particles randomly
      Particle P;
      P.sigma = sigma;
      for (int i=0; i<number; ++i) {
	P.position = size*vec2(drand48(), drand48());
	simData.addParticle(P);
      }
      // Set up sectorization
      Sectorization sectors(&simData);
      sectors.createVerletLists();
      
      // Look for minimum distance between particles
      RealType minDistance = sectors.getMinDistance();
      if (minDistance>sigma || minDistance<0) ++success;
    }

    // Record data
    data.push_back(pair<RealType, RealType>(phi, (RealType)success/trials));
  }

  cout << "data=" << data << ";\n";

  return 0;
}
