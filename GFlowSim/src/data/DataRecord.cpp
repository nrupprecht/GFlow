#include "DataRecord.hpp"
#include "../integrators/Integrator.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"
#include "../forces/ExternalForce.hpp"

namespace GFlow {
  DataRecord::DataRecord() : writeDirectory("RunData"), delay(1./15.), lastRecord(-delay), recIter(0), nsx(-1), nsy(-1), sdx(-1), sdy(-1), cutoff(-1), skinDepth(-1), recPos(false), recOption(1), recPerf(false), recMvRatio(false) {
    start_time = end_time = high_resolution_clock::now();
  };

  void DataRecord::startTiming() {
    start_time = high_resolution_clock::now();
  }

  void DataRecord::endTiming() {
    end_time = high_resolution_clock::now();
  }

  void DataRecord::markTime() {
    last_record = high_resolution_clock::now();
  }

  double DataRecord::getElapsedTime() const {
    return time_span(end_time, start_time);
  }

  void DataRecord::record(SimData* simData, RealType time) {
    // Return if not enough time has gone by
    if (time-lastRecord<delay) return;

    // Record performance
    if (recPerf) {
      high_resolution_clock::time_point current_time = high_resolution_clock::now();
      performanceRecord.push_back(pair<RealType,RealType>(time,time_span(current_time, last_record)));
      last_record = current_time;
    }

    // Record data
    lastRecord = time;

    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it = simData->getItPtr();
    int domain_size = simData->domain_size;

    // Record position data
    if (recPos) {
      vector<PData> positions;
      if (recOption==1) // Record pressure
	recordPressureData(simData, positions);
      else if (recOption==2) // Record by particle id
	recordByNumber(simData, positions);
      else if (recOption==3) // Record by number of verlet lists the particle is in
	recordByVerletList(simData, positions);
      else { // recOption==0 or default. Record nothing extra
	for (int i=0; i<domain_size; ++i)
	  if (it[i]>-1) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], 0));
      }
      positionRecord.push_back(positions);
    }

    // Record stat functions
    for (int i=0; i<statFunctions.size(); ++i) {
      StatFunc sf = statFunctions.at(i);
      statFunctionData.at(i).push_back(pair<RealType, RealType>(time, sf(simData)));
    }

    // Increment record counter
    recIter++;
  }

  void DataRecord::getSectorizationData(Sectorization* sectors) {
    if (sectors==nullptr) return;
    // Get width and height of sector grid
    nsx = sectors->nsx; 
    nsy = sectors->nsy;
    // Get sector's width and height
    sdx = sectors->sdx;
    sdy = sectors->sdy;
    // Get cutoff and skin depth
    cutoff = sectors->cutoff;
    skinDepth = sectors->skinDepth;
    // Get the number of neighbor lists
    auto verletList = sectors->verletList;
    numberOfVerletLists = verletList.size();
    // Get the average neighbors per verlet list
    RealType ave = 0;
    for (const auto& vl : verletList) ave += vl.size();
    if (numberOfVerletLists>0) ave /= static_cast<RealType>(numberOfVerletLists);
    else ave = -1;
    avePerVerletList = ave;
    
    // Get the number of occupied sectors
    auto sec = sectors->sectors;
    occupiedSectors = 0;
    avePerOccupiedSector = 0;
    for (int i=0; i<nsx*nsy; ++i) {
      if (sec[i].size()>0) {
	++occupiedSectors;
	avePerOccupiedSector += sec[i].size();
      }
    }
    avePerOccupiedSector /= (occupiedSectors>0 ? occupiedSectors : 1);
  }

  void DataRecord::writeData(SimData* simData) const {
    // Remove previously existing files
    system(("rm -rf "+writeDirectory).c_str());

    // Create the directory
    mkdir(writeDirectory.c_str(), 0777);

    // Create Positions directory
    if (recPos) {
      mkdir((writeDirectory+"/Pos").c_str(), 0777);
      // Write number
      printToCSV(writeDirectory+"/number.csv", vector<int>(1,recIter));
      // Write position data
      for (int i=0; i<positionRecord.size(); ++i)
	if (!printToCSV(writeDirectory+"/Pos/pos", positionRecord.at(i), i))
	  std::cerr << "Failed to print to [" << writeDirectory << "/Pos/pos" << i << "].\n";
      // Write wall data
      if (!printToCSV(writeDirectory+"/walls.csv", simData->getWalls()))
	std::cerr << "Failed to print walls to [" << writeDirectory << "/walls.csv].\n";
      // Write bounds data
      if (!printToCSV(writeDirectory+"/bnds.csv", vector<Bounds>(1,simData->getSimBounds())))
	std::cerr << "Failed to print bounds to [" << writeDirectory << "/bnds.csv].\n";
    }

    // Write stat function data to files
    if (!statFunctions.empty()) {
      mkdir((writeDirectory+"/StatData").c_str(), 0777);
      // Write the names of all the files that will be generated
      ofstream fout(writeDirectory+"/StatData/statNames.csv");
      if (fout.fail()) std::cerr << "Failed to open [" << writeDirectory << "/StatData/statNames.txt].\n";
      for (auto& name : statFunctionName) fout << name << "\n";
      fout.close();
      // Write stat data
      for (int i=0; i<statFunctionData.size(); ++i) {
	if (!printToCSV(writeDirectory+"/StatData/"+statFunctionName.at(i)+".csv", statFunctionData.at(i)))
	  std::cerr << "Failed to print to [" << writeDirectory << "/StatData/" << statFunctionName.at(i) << "].\n";
      }
    }
  }

  void DataRecord::writeRunSummary(SimData* simData, Integrator* integrator) const {
    std::ofstream fout(writeDirectory+"/run_summary.txt");
    if (fout.fail()) {
      // Write error message
      std::cerr << "Failed to open file [" << writeDirectory << "/run_summary.txt]." << endl;
      return;
    }
    // Print Header
    fout << "**********          SUMMARY          **********\n";
    fout << "**********  GFlow Granular Simulator **********\n";
    fout << "********** 2017, Nathaniel Rupprecht **********\n";
    fout << "***********************************************\n\n";
    // Print timing summary
    RealType elapsedTime = time_span(end_time, start_time);
    fout << "Timing and performance:\n";
    fout << "  - Simulated Time:           " << runTime << "\n";
    if (runTime!=actualTime) fout << "  - Actual time simulated:    " << actualTime;
    fout << "\n";
    fout << "  - Actual run Time:          " << elapsedTime << "\n";
    fout << "  - Ratio:                    " << runTime/elapsedTime << "\n";
    fout << "  - Inverse Ratio:            " << elapsedTime/runTime << "\n";
    fout << "\n";

    if (simData) {
      // Print simulation summary
      fout << "Simulation and space:\n";
      fout << "  - Dimensions:               " << simData->simBounds << "\n";
      fout << "  - Number of particles:      " << simData->domain_size << "\n";
      fout << "  - Ratio x Particles:        " << runTime/elapsedTime*simData->domain_size << "\n";
      fout << "\n";

      // Print forces
      if (!simData->externalForces.empty()) {
	fout << "External forces:\n";
	int i=1;
	for (const auto& f : simData->externalForces) {
	  fout << "  (" << i << ")\t" << f->summary() << "\n";
	  ++i;
	}
	fout << "\n";
      }

      // Print particle data
      if (0<simData->domain_size) writeParticleData(fout, simData);
    }

    if (integrator) {
      fout << "Integration:\n";
      fout << "  - Iterations:               " << integrator->iter << "\n";
      fout << "  - Time step (at end):       " << integrator->dt << "\n";
      auto vvint = reinterpret_cast<VelocityVerletIntegrator*>(integrator);
      if (vvint!=nullptr) {
	fout << "  - Average time step:        " << vvint->getAveTimeStep() << "\n";
	fout << "  - Average update delay:     " << vvint->getAveUpdateDelay() << "\n";
      }
      fout << "  - Time per iteration:       " << elapsedTime/integrator->iter << "\n";
      fout << "\n";
    }
    
    // Print the sectorization summary
    fout << "Sectorization summary (as of end of simulation):\n";
    fout << "  - Grid dimensions:          " << nsx << " x " << nsy << "\n";
    fout << "  - Grid lengths:             " << sdx << " x " << sdy << "\n";
    fout << "  - Number of verlet lists:   " << numberOfVerletLists << "\n";
    fout << "  - Average per verlet list:  " << (avePerVerletList>-1 ? toStr(avePerVerletList) : "---") << "\n";
    fout << "  - Cutoff:                   " << cutoff << "\n";
    fout << "  - Skin depth:               " << skinDepth << "\n";
    fout << "  - Occupied sectors:         " << occupiedSectors << "\n";
    fout << "  - Ave per occupied sector:  " << avePerOccupiedSector << "\n";
    // Close the stream
    fout.close();
  }

  void DataRecord::addStatFunction(StatFunc sf, string name) {
    // Add a place to store this function's data
    statFunctionData.push_back(vector<pair<RealType,RealType> >());

    // Add the stat function
    statFunctions.push_back(sf);

    // Add the stat function's name
    statFunctionName.push_back(name);
  }

  void DataRecord::push_mvRatio(RealType mv) {
    if (recMvRatio) movementRatioRecord.push_back(mv);
  }

  vector<pair<RealType, RealType> > DataRecord::getStatFunctionData(int index) const {
    if (index<0 || statFunctionData.size()<=index) throw BadStatFunction(index);
    return statFunctionData.at(index);
  }

  string DataRecord::getStatFunctionName(int index) const {
    if (index<0 || statFunctionData.size()<=index) throw BadStatFunction(index);
    return statFunctionName.at(index);
  }

  vector<pair<RealType, RealType> > DataRecord::getPerformanceRecord() const{
    return performanceRecord;
  }

  vector<RealType> DataRecord::getMoveRatioRecord() const {
    return movementRatioRecord;
  }
  
  void DataRecord::writeParticleData(std::ofstream& fout, SimData *simData) const {
    // We are here only if 0 < domain_size
    fout << "Particle Average Data:\n";

    // Find average sigma
    RealType sigma = 0;
    for (int i=0; i<simData->domain_size; ++i)
      if (-1<simData->it[i]) sigma += simData->sg[i];
    sigma /= simData->domain_size;
    fout << "  - Average sigma:            " << sigma << "\n";

    // Find average repulsion
    RealType repulsion = 0;
    for(int i=0; i<simData->domain_size; ++i)
      if (-1<simData->it[i]) repulsion += simData->rp[i];
    repulsion /= simData->domain_size;
    fout << "  - Average repulsion:        " << repulsion << "\n";

    // Find average dissipation
    RealType dissipation = 0;
    for(int i=0; i<simData->domain_size; ++i)
      if (-1<simData->it[i]) dissipation += simData->ds[i];
    dissipation /= simData->domain_size;
    fout << "  - Average dissipation:      " << dissipation << "\n";

    // Find average coefficient
    RealType coeff = 0;
    for(int i=0; i<simData->domain_size; ++i)
      if (-1<simData->it[i]) coeff += simData->cf[i];
    coeff /= simData->domain_size;
    fout << "  - Average coeff:            " << coeff << "\n";

    // Find average mass
    RealType mass = 0;
    for(int i=0; i<simData->domain_size; ++i)
      if (-1<simData->it[i]) mass += 1./simData->im[i];
    mass /= simData->domain_size;
    fout << "  - Average mass:             " << mass << "\n";

    // Find average density
    RealType density = 0;
    for(int i=0; i<simData->domain_size; ++i)
      if (-1<simData->it[i]) density += 1./simData->im[i]/(PI*sqr(simData->sg[i]));
    density /= simData->domain_size;
    fout << "  - Average density:          " << density << "\n";

    // Find average speed
    RealType speed = 0;
    for(int i=0; i<simData->domain_size; ++i)
      if (-1<simData->it[i]) speed += sqrt(sqr(simData->vx[i]) + sqr(simData->vy[i]));
    speed /= simData->domain_size;
    fout << "  - Average speed:            " << speed << "\n";
    
    // Find average omega
    RealType omega = 0;
    for(int i=0; i<simData->domain_size; ++i)
      if (-1<simData->it[i]) omega += simData->om[i];
    omega /= simData->domain_size;
    fout << "  - Average omega:            " << omega << "\n";

    // Section divide new line
    fout << "\n";
  }

  void DataRecord::recordPressureData(SimData* simData, vector<PData>& positions) const {
    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it = simData->getItPtr();
    int domain_size = simData->domain_size;

    // We will sort out particles with it<0 at the end
    vector<PData> pos;
    for (int i=0; i<domain_size; ++i) 
      pos.push_back(PData(px[i], py[i], sg[i], th[i], it[i], 0));
  
    // Get the force data
    ForceHandler *force = simData->getForceHandler();
    if (force && simData && simData->sectors) {
      // Get lists      
      const auto& verletList = simData->sectors->getVerletList();
      const auto& wallList   = simData->sectors->getWallList();
      // Get data
      force->pForcesRec(verletList, simData, pos);
      force->wForcesRec(wallList, simData, pos);
    }

    // Make sure we only have "real" particles (it>-1)
    for (int i=0; i<domain_size; ++i) 
      if (it[i]>-1) {
	PData pdata = pos.at(i);
	std::get<5>(pdata) /= 2*PI*sg[i]; // Convert to pressure
	positions.push_back(pdata);
      }
  }

  void DataRecord::recordByNumber(SimData* simData, vector<PData>& positions) const {
    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it = simData->getItPtr();
    int domain_size = simData->domain_size;

    for (int i=0; i<domain_size; ++i)
      if (it[i]>-1) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], i));
  }

  void DataRecord::recordByVerletList(SimData* simData, vector<PData>& positions) const {
    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it = simData->getItPtr();
    int domain_size = simData->domain_size;

    // We will sort out particles with it<0 at the end
    vector<PData> pos;
    for (int i=0; i<domain_size; ++i)
      pos.push_back(PData(px[i], py[i], sg[i], th[i], it[i], 0));
    
    if (simData && simData->sectors) {
      const auto& verletList = simData->sectors->getVerletList();
      for (auto vl : verletList)
	for (auto id : vl)
	  ++std::get<5>(pos.at(id));
    }
    
    // Make sure we only have "real" particles (it>-1)
    for (int i=0; i<domain_size; ++i)
      if (it[i]>-1) positions.push_back(pos.at(i));
  }

}
