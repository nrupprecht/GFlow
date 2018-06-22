/*
 * Author: Nathaniel Rupprecht
 * Start Data: June 15, 2018
 *
 */

#include "../../include/ArgParse.hpp"
#include "../control/Sectorization.hpp"
#include "../creation/FileParser.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"
#include "../data/DataRecord.hpp"

using namespace GFlow;

int main(int argc, char** argv) {
  // Variables
  int trials = 1000;
  RealType minPhi = 0;
  RealType maxPhi = 0.1;
  RealType numPhi = 10;
  RealType sigma = 0.01;
  RealType size = 1.;

  // Set up some pointers to objects we'll need
  Integrator *integrator = nullptr;
  DataRecord *dataRecord = nullptr;
  SimData *simData = nullptr;
  
  // Read some simulation options
  ArgParse parser(argc, argv);
  parser.get("trials", trials);
  parser.get("minPhi", minPhi);
  parser.get("maxPhi", maxPhi);
  parser.get("numPhi", numPhi);
  parser.get("size", size);

  // Set up file parser
  FileParser fileParser;
  fileParser.set(argc, argv);
  fileParser.setArgParse(&parser);
  fileParser.setDataRecord(dataRecord);

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

    // Set up particles to be non-overlapping
    fileParser.parse("samples/square.cfg", simData, integrator);

    int success = 0;
    for (int tr = 0; tr < trials; ++tr) {
      // Set up simulation data
      SimData simData(bounds, bounds);
      simData.setWrapX(true);
      simData.setWrapY(true);
      
      // Add particles randomly
      Particle P;
      P.sigma = sigma;
      P.position = size*vec2(drand48(), drand48());

      // Set up sectorization
      Sectorization sectors(&simData);
      sectors.createVerletLists();

      // Let particles relax

      // Insert a new particle

      // Check if it overlaps with the others
    }

    // Record data
    data.push_back(pair<RealType, RealType>(phi, (RealType)success/trials));
  }

  cout << "data=" << data << ";\n";

  return 0;
}
