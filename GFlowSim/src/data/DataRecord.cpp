#include "DataRecord.hpp"
#include "../integrators/Integrator.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"
#include "../forces/ExternalForce.hpp"

namespace GFlow {
  DataRecord::DataRecord() : doRecord(true), simData(nullptr), sectors(nullptr), forceHandler(nullptr), writeDirectory("RunData"), delay(1./15.), lastRecord(-2.*delay), animationDelay(1./15.), lastAnimate(-2.*animationDelay), recIter(0), nsx(-1), nsy(-1), sdx(-1), sdy(-1), cutoff(-1), skinDepth(-1), recPos(false), lowerSizeLimit(0), recOption(1), stripeHeight(3), recBulk(false), recBulkOutline(false), recDisplacementField(false), displacementSnapshot(false), trackDisplacement(false), recPressField(false), recPerf(false), statRecPerf(-1), recMvRatio(false), statRecMvRatio(-1), recDt(false), statRecDt(-1), recDelay(false), statRecDelay(-1), center(false) {
    start_time = end_time = high_resolution_clock::now();
  };

  void DataRecord::initialize() {
    int end = statFunctions.size();
    // Put in other statistics (not from normal stat functions) after the normal stat functions
    if (recPerf) {
      statFunctionData.push_back(vector<RPair>());
      statFunctionName.push_back("perf");
      statRecPerf = end;
      ++end;
    }
    if (recMvRatio) {
      statFunctionData.push_back(vector<RPair>());
      statFunctionName.push_back("mvRatio");
      statRecMvRatio = end;
      ++end;
    }
    if (recDt) {
      statFunctionData.push_back(vector<RPair>());
      statFunctionName.push_back("dt");
      statRecDt = end;
      ++end;
    }
    if (recDelay) {
      statFunctionData.push_back(vector<RPair>());
      statFunctionName.push_back("delay");
      statRecDelay = end;
      ++end;
    }
    if (recDistance) {
      statFunctionData.push_back(vector<RPair>());
      statFunctionName.push_back("distance");
      statRecDistance = end;
      ++end;
    }
  }

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

  RealType DataRecord::getRatio() const {
    return actualTime/getElapsedTime();
  }

  void DataRecord::record(RealType time) {
    // Checks
    if (!doRecord) return;

    // Record position data - on a different timer than stat data things
    if (recPos && time-lastAnimate>animationDelay) {
      vector<PData> positions;
      getPositionData(positions, recOption);
      positionRecord.push_back(positions);
      lastAnimate = time;
    }

    // Return if not enough time has gone by
    if (time-lastRecord<delay) return;
    // Record performance
    if (recPerf) {
      high_resolution_clock::time_point current_time = high_resolution_clock::now();
      statFunctionData.at(statRecPerf).push_back(RPair(time,time_span(current_time, last_record)));
      last_record = current_time;
    }

    // Record data
    lastRecord = time;

    // Initialize displacement tracking, if neccessary
    if (recIter==0) {
      // Set our own initial position record (unneccessary now?)
      RealType *px = simData->getPxPtr();
      RealType *py = simData->getPyPtr();
      int domain_end = simData->domain_end;    
      initialPositions = vector<vec2>(domain_end);
      // Record
      for (int i=0; i<domain_end; ++i) initialPositions.at(i) = vec2(px[i], py[i]);
    }

    // Record position data
    if (recPos && time-lastAnimate>animationDelay) {
      vector<PData> positions;
      getPositionData(positions, recOption);
      positionRecord.push_back(positions);
      lastAnimate = time;
    }

    // Record bulk data
    if (recBulk || recBulkOutline) {
      Bounds region;
      if (center) {
	RealType max_posx =StatFunc_MaxR_PosX(simData);
	RealType max_posy = StatFunc_MaxR_PosY(simData);
	RealType maxSigma = StatFunc_MaxR_Sigma(simData);
	region = Bounds(max_posx-3*maxSigma, max_posx+3*maxSigma, max_posy-maxSigma, max_posy+8*maxSigma);
      }
      else region = simData->getSimBounds();
      vector<RealType> volumes;
      vector<pair<vec2,vec2> > lines;
      getBulkData(region, volumes, bulkField, lines);
      if (recBulkOutline) bulkOutline.push_back(lines);
      if (recBulk)        bulkVolumes.push_back(volumes);
    }

    // Record pressure data 
    if (recPressField) {
      Bounds region;
      if (center) {
	RealType max_posx = StatFunc_MaxR_PosX(simData);
	RealType max_posy = StatFunc_MaxR_PosY(simData);
	RealType maxSigma = StatFunc_MaxR_Sigma(simData);
	region = Bounds(max_posx-6*maxSigma, max_posx+6*maxSigma, max_posy-6*maxSigma, max_posy+18*maxSigma);
      }
      else region = simData->getSimBounds();
      getPressureData(region, pressField);
    }

    // Record vortex data
    if (recVortex) {
      Bounds region;
      if (center) {
        RealType max_posx = StatFunc_MaxR_PosX(simData);
        RealType max_posy = StatFunc_MaxR_PosY(simData);
        RealType maxSigma = StatFunc_MaxR_Sigma(simData);
        region = Bounds(max_posx-6*maxSigma, max_posx+6*maxSigma, max_posy-6*maxSigma, max_posy+18*maxSigma);
      }
      else region = simData->getSimBounds();
      ScalarField vortex;
      getVortexData(region, vortex);
      vortexRecord.push_back(vortex);
    }

    // Record displacement field
    if (recDisplacementField) displacementFieldRecord.push_back(createDisplacementField(true));

    // Record stat functions
    for (int i=0; i<statFunctions.size(); ++i) {
      StatFunc sf = statFunctions.at(i);
      statFunctionData.at(i).push_back(pair<RealType, RealType>(time, sf(simData)));
    }

    // Record stat plots
    for (int i=0; i<statPlots.size(); ++i) {
      StatPlot sp = statPlots.at(i);
      sp(simData, statPlotData.at(i), statPlotBounds.at(i));
    }

    // Record average distance
    if (recDistance && statRecDistance>-1) {
      // Find the average distance
      RealType R = 0;
      RealType *px = simData->getPxPtr();
      RealType *py = simData->getPyPtr();
      int *it = simData->getItPtr();
      int domain_end = simData->getDomainEnd(), domain_size = simData->getDomainSize();
      if (domain_size>0) {
	for (int i=0; i<domain_end; ++i) {
	  if (-1<it[i]) R += sqrt(sqr(simData->getDisplacement(vec2(px[i], py[i]), initialPositions.at(i))));
	}
	R /= domain_size;
	statFunctionData.at(statRecDistance).push_back(RPair(time, R));
      }
    }

    // Increment record counter
    recIter++;
  }

  void DataRecord::record(VelocityVerletIntegrator* integrator, RealType time) {
    // Record move ratio
    if (recMvRatio && statRecMvRatio>-1) statFunctionData.at(statRecMvRatio).push_back(RPair(time, integrator->getMvRatio()));

    // Record time step
    if (recDt && statRecDt>-1) statFunctionData.at(statRecDt).push_back(RPair(time, integrator->getDt()));

    // Record update delay
    if (recDelay && statRecDelay>-1) statFunctionData.at(statRecDelay).push_back(RPair(time, integrator->getUpdateDelay()));
 }

  void DataRecord::getSectorizationData() {
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

  void DataRecord::writeData() const {
    // Remove previously existing files
    system(("rm -rf "+writeDirectory).c_str());

    // Create the directory
    mkdir(writeDirectory.c_str(), 0777);

    // Create Positions directory
    if (recPos) {
      // Print position data
      printToDirectory(writeDirectory+"/Pos", "pos", positionRecord);

      // Write wall data
      if (!printToCSV(writeDirectory+"/walls.csv", simData->getWalls()))
	std::cerr << "Failed to print walls to [" << writeDirectory << "/walls.csv].\n";
    }
     
    // Write bounds data
    if (!printToCSV(writeDirectory+"/bnds.csv", vector<Bounds>(1,simData->getSimBounds())))
      std::cerr << "Failed to print bounds to [" << writeDirectory << "/bnds.csv].\n";

    // Write stat function data to files
    if (!statFunctionData.empty()) {
      mkdir((writeDirectory+"/StatData").c_str(), 0777);
      // Write the names of all the files that will be generated
      ofstream fout(writeDirectory+"/StatData/statNames.csv");
      if (fout.fail()) std::cerr << "Failed to open [" << writeDirectory << "/StatData/statNames.txt].\n";
      for (auto& name : statFunctionName) fout << name << "\n";
      fout.close();
      // Write stat function data
      for (int i=0; i<statFunctionData.size(); ++i) {
	if (!printToCSV(writeDirectory+"/StatData/"+statFunctionName.at(i)+".csv", statFunctionData.at(i)))
	  std::cerr << "Failed to print to [" << writeDirectory << "/StatData/" << statFunctionName.at(i) << ".csv].\n";
      }
    }

    // Write stat plot data to files
    if (!statPlotData.empty()) {
      mkdir((writeDirectory+"/StatPlotData").c_str(), 0777);
      // Write the names of all the files that will be generated
      ofstream fout(writeDirectory+"/StatPlotData/statNames.csv");
      if (fout.fail()) std::cerr << "Failed to open [" << writeDirectory << "/StatPlotData/statNames.txt].\n";
      for (auto& name : statPlotName) fout << name << "\n";
      fout.close();
      // Write stat plot data
      for (int i=0; i<statPlotData.size(); ++i) {
	if (!printToCSV(writeDirectory+"/StatPlotData/"+statPlotName.at(i)+".csv", statPlotData.at(i)))
	  std::cerr << "Failed to print to [" << writeDirectory << "/StatPlotData/" << statPlotName.at(i) << ".csv].\n";
      }
    }

    // Write displacement data
    if (trackDisplacement) writeDisplacementData();
    if (displacementSnapshot) writeDisplacementField();
    if (recDisplacementField) {
      // Create a directory for the field frames
      mkdir((writeDirectory+"/DisplacementField").c_str(), 0777);
      // Write the fields
      for (int i=0; i<displacementFieldRecord.size(); ++i) {
	displacementFieldRecord.at(i).printToCSV(writeDirectory+"/DisplacementField/dispField"+toStr(i)+".csv");
      }
    }
    if (simData && recDisplacementColumns) {
      auto sb = simData->getSimBounds();
      int bins = (sb.right - sb.left)/0.2;
      writeColumnDisplacement(bins);
    }

    // Write bulk data
    if (recBulk) {
      bulkField.printToCSV(writeDirectory+"/bulkField.csv");
      printToCSV(writeDirectory+"/bulkVolumes.csv", bulkVolumes);
    }

    // Writebulk outlines
    if (recBulkOutline) printToDirectory(writeDirectory+"/bulkOutline", "bulk", bulkOutline);

    // Write pressure data
    if (recPressField) pressField.printToCSV(writeDirectory+"/pressField.csv", 1./static_cast<RealType>(recIter));

    // Write vortex data
    if (recVortex) {
      // Create vortex directors
      mkdir((writeDirectory+"/Vortex").c_str(), 0777);
      // Print frames
      int index = 0;
      for (const auto &f : vortexRecord) {
	f.printToCSV(writeDirectory+"/Vortex/vortex"+toStr(index)+".csv");
	++index;
      }
    }
  }

  void DataRecord::writeRunSummary(Integrator* integrator) const {
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
    // Print command
    if (!command.empty()) {
      fout << "Command:\n  ";
      for (auto& c : command) fout << c << " ";
      fout << "\n\n";
    }
    // Print timing summary
    RealType elapsedTime = time_span(end_time, start_time);
    fout << "Timing and performance:\n";
    fout << "  - Set up time:              " << setupTime << "\n";
    fout << "  - Time simulated:           " << actualTime << "\n";
    fout << "  - Requested Time:           " << runTime << "\n";
    fout << "  - Run Time:                 " << elapsedTime;
    if (elapsedTime>60) fout << " ( h:m:s - " <<printAsTime(elapsedTime) << " )";
    fout << "\n";
    fout << "  - Ratio:                    " << actualTime/elapsedTime << "\n";
    fout << "  - Inverse Ratio:            " << elapsedTime/actualTime << "\n";
    fout << "\n";

    if (simData) {
      // Print simulation summary
      fout << "Simulation and space:\n";
      fout << "  - Dimensions:               " << simData->simBounds << "\n";
      fout << "  - Wrapping:                 " << simData->getWrapX() << ", " << simData->getWrapY() << "\n";
      fout << "  - Number of particles:      " << simData->domain_size << "\n";
      fout << "  - Ratio x Particles:        " << runTime/elapsedTime*simData->domain_size << "\n";
      RealType vol = 0;
      auto particles = simData->getParticles();
      for (const auto& p : particles) vol += sqr(p.sigma);
      vol *= PI;
      RealType phi = vol/simData->getSimBounds().volume();
      fout << "  - Packing fraction:         " << phi << "\n";
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
      if (0<simData->domain_size) writeParticleData(fout);

      // Print checks
      if (0<simData->domain_size) writeParticleChecks(fout);
    }

    if (integrator) {
      fout << "Integration:\n";
      fout << "  - Iterations:               " << integrator->iter << "\n";
      fout << "  - Time per iteration:       " << elapsedTime / integrator->iter << "\n";
      fout << "  - Time step (at end):       " << integrator->dt << "\n";
      auto vvint = dynamic_cast<VelocityVerletIntegrator*>(integrator);
      if (vvint!=nullptr) {
	fout << "  - Average time step:        " << vvint->getAveTimeStep() << "\n";
	fout << "  - Average update delay:     " << vvint->getAveUpdateDelay() << "\n";
      }
      fout << "  - Time per iteration:       " << elapsedTime/integrator->iter << "\n";
      fout << "\n";
    }
    
    // Print the sectorization summary
    if (sectors) {
      fout << "Sectorization summary (as of end of simulation):\n";
      int nsx = sectors->getNSX(), nsy = sectors->getNSY();
      
      fout << "  - Grid dimensions:          " << nsx << " x " << nsy << "\n";
      fout << "  - Total sectors:            " << nsx * nsy << "\n";
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
  }

  void DataRecord::setCommand(int argc, char** argv) {
    command.clear();
    for (int i=0; i<argc; ++i) command.push_back(argv[i]);
  }

  void DataRecord::resetTimers() {
    lastRecord = -2.*delay;
    lastAnimate = -2.*animationDelay;
  }

  void DataRecord::addStatFunction(StatFunc sf, string name) {
    // Add the stat function
    statFunctions.push_back(sf);

    // Add a place to store this function's data
    statFunctionData.push_back(vector<pair<RealType,RealType> >());

    // Add the stat function's name
    statFunctionName.push_back(name);
  }

  void DataRecord::addStatPlot(StatPlot sp, RPair bounds, int bins, string name) {
    // Add the stat plot
    statPlots.push_back(sp);

    // Add a place to store this plot's data
    vector<RPair> plot(bins);
    RealType dr = (bounds.second-bounds.first)/static_cast<RealType>(bins);
    for (int i=0; i<bins; ++i) {
      plot.at(i).first  = bounds.first + i*dr;
      plot.at(i).second = 0;
    }
    statPlotData.push_back(plot);

    // Add the stat plot's bounds
    statPlotBounds.push_back(bounds);

    // Add the stat plot's name
    statPlotName.push_back(name);
  }

  void DataRecord::getPositionData(vector<PData>& positions, int option) {
    if (option==1) // Record pressure
      recordByPressure(positions);
    else if (option==2) // Record by particle id
      recordByNumber(positions);
    else if (option==3) // Record by number of verlet lists the particle is in
      recordByVerletList(positions);
    else if (option==4) // Record by velocity
      recordByVelocity(positions);
    else if (option==5) // Record by displacement
      recordByDisplacement(positions);
    else if (option==6) // Record by height, in colored stripes
      recordByHeight(positions);
    else if (option==7) // Color marked particles
      recordByMarked(positions);
    else { // recOption==0 or default. Record nothing extra
      RealType *px = simData->getPxPtr();
      RealType *py = simData->getPyPtr();
      RealType *sg = simData->getSgPtr();
      RealType *th = simData->getThPtr();
      int *it = simData->getItPtr();
      int domain_end = simData->domain_end;
      for (int i=0; i<domain_end; ++i)
	if (it[i]>-1 && lowerSizeLimit<sg[i]) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], 0));
    }
  }

  vector<pair<RealType, RealType> > DataRecord::getStatFunctionData(int index) const {
    if (index<0 || statFunctionData.size()<=index) throw BadStatFunction(index);
    return statFunctionData.at(index);
  }

  string DataRecord::getStatFunctionName(int index) const {
    if (index<0 || statFunctionData.size()<=index) throw BadStatFunction(index);
    return statFunctionName.at(index);
  }

  void DataRecord::writeParticleData(std::ofstream& fout) const {
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

  void DataRecord::writeParticleChecks(std::ofstream& fout) const {
    // We are here only if 0 < domain_size
    fout << "Particle Checks:\n";
    
    // Find max and min particle positions
    RealType minX(1e9), maxX(-1e9), minY(1e9), maxY(-1e9);
    for (int i=0; i<simData->domain_size; ++i) {
      if (simData->it[i]>-1) {
	if (maxX < simData->px[i]) maxX = simData->px[i];
	if (simData->px[i]<minX)   minX = simData->px[i];
	if (maxY < simData->py[i]) maxY = simData->py[i];
	if (simData->py[i]<minY)   minY = simData->py[i];
      }
    }
    fout << "  - (Min x, Max x):    " << minX << ", " << maxX << "\n";
    fout << "  - (Min y, Max y):    " << minY << ", " << maxY << "\n";
    fout << "\n";
  }

  void DataRecord::recordByPressure(vector<PData>& positions) const {
    // Checks
    if (sectors==nullptr || forceHandler==nullptr) return;

    // Get arrays
    auto px = simData->getPxPtr();
    auto py = simData->getPyPtr();
    auto sg = simData->getSgPtr();
    auto th = simData->getThPtr();
    auto it = simData->getItPtr();
    int domain_end = simData->getDomainEnd();

    // Add data to the vector, except for special data (last element)
    vector<PData> pos;
    for (int i=0; i<domain_end; ++i)
      pos.push_back(PData(px[i], py[i], sg[i], th[i], it[i], 0));

    // Get lists
    const auto& verletList = sectors->getVerletList();
    const auto& wallList   = sectors->getWallList();

    // Get data
    forceHandler->pForcesRec(verletList, simData, pos);
    forceHandler->wForcesRec(wallList, simData, pos);

    // Make sure we only have "real" particles (it>-1)
    for (int i=0; i<domain_end; ++i)
      if (it[i]>-1 && lowerSizeLimit<sg[i]) {
        PData pdata = pos.at(i);
	std::get<5>(pdata) /= 2*PI*sg[i]; // Convert to pressure
        positions.push_back(pdata);
      }
  }

  void DataRecord::recordByNumber(vector<PData>& positions) const {
    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it = simData->getItPtr();
    int domain_size = simData->domain_size;

    for (int i=0; i<domain_size; ++i)
      if (it[i]>-1 && lowerSizeLimit<sg[i]) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], i));
  }

  void DataRecord::recordByVerletList(vector<PData>& positions) const {
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
    
    if (simData && sectors) {
      const auto& verletList = sectors->getVerletList();
      for (auto vl : verletList)
	for (auto id : vl)
	  ++std::get<5>(pos.at(id));
    }
    
    // Make sure we only have "real" particles (it>-1)
    for (int i=0; i<domain_size; ++i)
      if (it[i]>-1 && lowerSizeLimit<sg[i]) positions.push_back(pos.at(i));
  }

  void DataRecord::recordByVelocity(vector<PData>& positions) const {
    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    int *it = simData->getItPtr();
    int domain_end = simData->domain_end;

    // We will sort out particles with it<0 at the end
    for (int i=0; i<domain_end; ++i)
      if (it[i]>-1 && lowerSizeLimit<sg[i]) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], sqrt(sqr(vx[i])+sqr(vy[i]))));
  }

  void DataRecord::recordByDisplacement(vector<PData>& positions) const {
    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    int *it = simData->getItPtr();
    int domain_end = simData->domain_end;

    // We will sort out particles with it<0 at the end
    for (int i=0; i<domain_end; ++i)
      if (it[i]>-1) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], sqrt(sqr(initialPositions.at(i)-vec2(px[i],py[i])))));
  }

  void DataRecord::recordByHeight(vector<PData>& positions) const {
    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    int *it = simData->getItPtr();
    int domain_end = simData->domain_end;

    // We will sort out particles with it<0 at the end
    for (int i=0; i<domain_end; ++i) {
      int col = static_cast<int>(initialPositions.at(i).y / stripeHeight) % 2;
      if (it[i]>-1) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], col));
    }
  }

  void DataRecord::recordByMarked(vector<PData>& positions) const {
    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    int *it = simData->getItPtr();
    int domain_end = simData->domain_end;

    // We will sort out particles with it<0 at the end
    for (int i=0; i<domain_end; ++i) {
      // Unordered set returns an iterator to the element if it is found
      RealType col = marked.find(i)!=marked.end() ? 1. : 0.;
      // Push back the data
      if (it[i]>-1) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], col));
    }
  }

  void DataRecord::getBulkData(const Bounds& region, vector<RealType>& volumes, ScalarField& field, vector<pair<vec2,vec2> >& lines, RealType resolution, RealType dr, RealType lowerVCut, RealType upperVCut) const {
    int nsx = (region.right - region.left)/resolution, nsy = (region.top - region.bottom)/resolution;
    RealType sx = (region.right - region.left)/nsx, sy = (region.top - region.bottom)/nsy;
    ++nsx;
    ++nsy;
    list<int> *bulkSectors = new list<int>[nsx*nsy];
    // Fill sectors
    int i=0;
    auto allParticles = simData->getParticles();
    for (const auto &p : allParticles) {
      int sec_x = (p.position.x - region.left)/sx;
      int sec_y = (p.position.y - region.bottom)/sy;
      if (-1<sec_x && sec_x<nsx && -1<sec_y && sec_y<nsy)
	bulkSectors[nsx*sec_y+sec_x].push_back(i);
      // If part of the particle sticks into the region, place the particle in the closest sector. This ignores particles "diagonal" to the region that stick into the region
      else if (region.left<p.position.x+p.sigma   && -1<sec_y && sec_y<nsy)
	bulkSectors[nsx*sec_y+0].push_back(i); // Left
      else if (p.position.x-p.sigma<region.right  && -1<sec_y && sec_y<nsy)
	bulkSectors[nsx*sec_y+(nsx-1)].push_back(i); // Right
      else if (region.bottom<p.position.y+p.sigma && -1<sec_x && sec_x<nsx)
	bulkSectors[sec_x].push_back(i); // Bottom
      else if (p.position.y-p.sigma<region.top    && -1<sec_x && sec_x<nsx)
	bulkSectors[nsx*(nsy-1)+sec_x].push_back(i); // Top
      ++i;
    }
    // Find the particle of maximum radius
    double maxR = 0;
    for (const auto &p : allParticles)
      if (p.sigma>maxR) maxR = p.sigma;
    // Check which sector centers are covered by particles
    int *array = new int[nsx*nsy];
    int sweepX = (maxR+dr)/sx, sweepY = (maxR+dr)/sy;
    for (int y=0; y<nsy; ++y)
      for (int x=0; x<nsx; ++x) {
	// Check whether center of sector is covered by any particles
	array[nsx*y+x] = nsx*y+x;
	bool done = false;
	vec2 pos((x+0.5)*sx+region.left, (y+0.5)*sy+region.bottom);
	// If outside the simulation bounds, there are no bubbles here
	if (!simData->simBounds.contains(pos)) {
	  array[nsx*y+x] = -1;
	  continue;
	}
	// Sweep through sectors
	int startX = max(0, x-sweepX), endX = min(nsx-1, x+sweepX);
	int startY = max(0, y-sweepY), endY = min(nsy-1, y+sweepY);
	// Check your own sector first for speed's sake
	for (const auto index : bulkSectors[nsx*y+x]) {
	  vec2 dR = pos-allParticles.at(index).position;
	  RealType R = dr+allParticles.at(index).sigma;
	  if (sqr(dR)<sqr(R)) {
	    done = true;
	    array[nsx*y+x] = -1;
	    break;
	  }
	}
	// Search the other sectors
	for (int Y=startY; Y<endY && !done; ++Y) {
	  for (int X=startX; X<endX && !done; ++X) {
	    if (X!=x || Y!=y)
	      for (const auto index : bulkSectors[nsx*Y+X]) {
		vec2 dR= pos-allParticles.at(index).position;
		RealType R = dr+allParticles.at(index).sigma;
		if (sqr(dR)<sqr(R)) {
		  done = true;
		  array[nsx*y+x] = -1;
		  break;
		}
	      }
	    if (done) break;
	  }
	  if (done) break;
	}
      }
    // Unite bubbles
    unite(array, nsx, nsy);
    // Find all bubbles that border an edge
    std::set<int> edgeBubbles;
    auto contains = [&] (int i) { return edgeBubbles.find(i)!=edgeBubbles.end(); };
    auto add = [&] (int i) {
      int h = getHead(array, i);
      if (h<0 || contains(h));
      else edgeBubbles.insert(h);
    };
    for (int x=0; x<nsx; ++x) add(x); // Bottom
    for (int x=0; x<nsx; ++x) add(nsx*nsy-x-1); // Top
    for (int y=0; y<nsy; ++y) add(nsx*y); // Left
    for (int y=0; y<nsy; ++y) add(nsx*y+nsx-1); // Right
    // Point all sectors to their head
    for (int k=0; k<nsx*nsy; ++k) array[k] = getHead(array,k);
    // "Erase" bubbles that touch the bounds
    for (int k=0; k<nsx*nsy; ++k)
      if (contains(array[k])) array[k] = -1;
    // Collect all head nodes
    std::map<int, int> labels; // Maps old label to new label
    std::map<int, int> volCount; // Maps (new) label to number of cells in the bubble
    int lab = 0;
    for (int k=0; k<nsx*nsy; ++k)
      if (array[k]!=-1)
	// If we have not already recorded this head
	if (labels.find(array[k])==labels.end()) {
	  labels.insert(pair<int, int>(array[k], lab));
	  volCount.insert(pair<int, int>(lab, 0));
	  ++lab;
	}
    // If there are no empty volumes
    if (labels.empty()) {
      delete [] bulkSectors;
      delete [] array;
      return;
    }
    // Count volumes
    for (int k=0; k<nsx*nsy; ++k) {
      if (array[k]!=-1) { // Search for the proper index
	int j=0;
	// Increment volume counter
	auto it = volCount.find(array[k]);
	if (it!=volCount.end())
	  ++it->second;
      }
    }
    // Record volume sizes
    i=0;
    volumes.clear();
    for (const auto c : volCount) {
      double vol = sx*sy*c.second;
      if (lowerVCut<vol && vol<upperVCut) volumes.push_back(vol);
    }
    // Remove volumes that are to small
    for (int y=0; y<nsy; ++y)
      for (int x=0; x<nsx; ++x) {
	int j = array[nsx*y+x];
	if (-1<j) {
	  auto it = volCount.find(j);
	  if (sx*sy*it->second<upperVCut) array[nsx*y+x] = -1;
	}
      }
    // Record data in scalar field
    if (!volumes.empty()) {
      // Set up field
      if (field.empty()) {
	field.setBounds(region);
	field.setResolution(sx);
	field.setPrintPoints(200);
      }
      // Update
      try {
	for (int y=0; y<field.getNSY(); ++y)
	  for (int x=0; x<field.getNSX(); ++x)
	    field.at(x,y) += (array[nsx*y+x]>-1 ? 1 : 0);
      }
      catch (ScalarField::FieldOutOfBounds bnds) {
	cout << "Current dims: " << nsx << ", " << nsy << endl;
	cout << "Tried to access " << bnds.x << ", " << bnds.y << endl;
	cout << "Field dims: " << field.getNSX() << ", " << field.getNSY() << endl;
      }
    }

    // Create outline
    if (recBulkOutline) createOutline(array, nsx, nsy, sx, sy, region, lines);
    
    // Sort sizes
    std::sort(volumes.begin(), volumes.end());
    // Clean up and return
    delete [] bulkSectors;
    delete [] array;
  }

  inline void DataRecord::unite(int *array, int nsx, int nsy) const {
    // Find out which volumes are connected
    for (int y=0; y<nsy; ++y)
      for (int x=0; x<nsx; ++x) {
	// Check the sectors around you, point yourself to the largest head
	int h0=-1, h1=-1, h2=-1, h3=-1;
	if (array[nsx*y+x]!=-1) {
	  int head = getHead(array, nsx*y+x);
	  if (0<x && -1<array[nsx*y+x-1]) { // Left
	    h0 = getHead(array, nsx*y+(x-1));
	    if (h0<head) head = h0;
	  }
	  if (x+1<nsx && -1<array[nsx*y+x+1]) { // Right
	    h1 = getHead(array, nsx*y+(x+1));
	    if (h1<head) head = h1;
	  }
	  if (0<y && -1<array[nsx*(y-1)+x]) { // Bottom
	    h2 = getHead(array, nsx*(y-1)+x);
	    if (h2<head) head = h2;
	  }
	  if (y+1<nsy && -1<array[nsx*(y+1)+x]) { // Top
	    h3 = getHead(array, nsx*(y+1)+x);
	    if (h3<head) head = h3;
	  }
	  // Set your head and the heads of all your neighbors as well
	  array[nsx*y+x] = head;
	  if (-1<h0) array[h0] = head;
	  if (-1<h1) array[h1] = head;
	  if (-1<h2) array[h2] = head;
	  if (-1<h3) array[h3] = head;
	}
      }
  }

  inline void DataRecord::createOutline(int *array, int nsx, int nsy, RealType sx, RealType sy, Bounds region, vector<pair<vec2,vec2> >& lines) const {
    // Create a lambda function to edge detect
    auto isEdge = [&] (int x, int y) {
      if (x<0 || nsx<=x || y<0 || nsy<=y || array[nsx*y+x]!=-1) return false;
      if (x+1<nsx && array[nsx*y+x+1]>-1)  return true;
      else if (0<x && array[nsx*y+x-1]>-1) return true;
      else if (y+1<nsy && array[nsx*(y+1)+x]>-1) return true;
      else if (0<y && array[nsx*(y-1)+x]>-1) return true;
      return false;
    };

    typedef pair<vec2,vec2> VPair;
    
    lines.clear();
    for (int y=0; y<nsy; ++y)
      for (int x=0; x<nsx; ++x) {
	if (isEdge(x,y)) {
	  vec2 pos((x+0.5)*sx+region.left, (y+0.5)*sy+region.bottom);
	  if (isEdge(x+1,y-1)) lines.push_back(VPair(pos, pos+vec2(sx,-sy)));
	  if (isEdge(x+1,y))   lines.push_back(VPair(pos, pos+vec2(sx,0)));
	  if (isEdge(x+1,y+1)) lines.push_back(VPair(pos, pos-vec2(sx,sy)));
	  if (isEdge(x,y+1))   lines.push_back(VPair(pos, pos+vec2(0,sy)));
	}
      }
  }

  inline int DataRecord::getHead(int* array, int index) const {
    if (index<0) return -1;
    while (array[index]!=index && -1<array[index]) index = array[index];
    return index;
  }

  inline void DataRecord::writeDisplacementData() const {
    // Get data
    int domain_end = simData->getDomainEnd();
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it      = simData->getItPtr();
    // We omit the largest particle
    RealType maxSg = StatFunc_MaxR_Sigma(simData);
    // Collect data
    vector<PData> displacementI;
    vector<PData> displacementF;
    for (int i=0; i<domain_end; ++i) {
      RealType disp = sqrt(sqr(initialPositions.at(i) - vec2(px[i], py[i])));
      if (sg[i]!=maxSg) { // ALWAYS EXCLUDES THE LARGEST PARTICLE
	PData dataI(initialPositions.at(i).x, initialPositions.at(i).y, sg[i], th[i], it[i], disp);
	PData dataF(px[i], py[i], sg[i], th[i], it[i], disp);
	displacementI.push_back(dataI);
	displacementF.push_back(dataF);
      }
    }
    // Write displacement to file
    if (!printToCSV(writeDirectory+"/displacementI.csv", displacementI))
      std::cerr << "Failed to print to [" << writeDirectory << "/displacementI.csv].\n";
    if (!printToCSV(writeDirectory+"/displacementF.csv", displacementF)) 
      std::cerr << "Failed to print to [" << writeDirectory << "/displacementF.csv].\n";
  }

  inline ScalarField DataRecord::createDisplacementField(bool finalPos) const {
    // Get data
    int domain_end = simData->getDomainEnd();
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it      = simData->getItPtr();
    // We may omit the largest particle
    RealType maxSg = StatFunc_MaxR_Sigma(simData);
    // Create field
    ScalarField field;
    Bounds bounds = simData->getSimBounds();
    field.setBounds(bounds);
    field.setResolution(0.1);
    field.setPrintPoints(400);
    RealType idx = 1./field.getDx(), idy = 1./field.getDy();
    int nsx = field.getNSX(), nsy = field.getNSY();
    int *count = new int[nsx*nsy];
    memset(count, 0, nsx*nsy*sizeof(int)); // Zero count array
    // Collect data
    for (int i=0; i<domain_end; ++i) {
      if (sg[i]==maxSg) continue;
      RealType disp = sqrt(sqr(simData->getDisplacement(initialPositions.at(i), vec2(px[i], py[i]))));
      // Bin data
      int sec_x = (initialPositions.at(i).x-bounds.left)*idx;
      int sec_y = (initialPositions.at(i).y-bounds.bottom)*idy;
      ++count[sec_y*nsx+sec_x];
      field.at(sec_x, sec_y) += disp;
    }
    for (int y=0; y<nsy; ++y)
      for (int x=0; x<nsx; ++x)
	if (count[y*nsx+x]>0) field.at(x,y) /= (RealType)count[y*nsx+x];
    // Return field
    return field;
  }

  inline void DataRecord::writeDisplacementField() const {
    ScalarField field = createDisplacementField(true);
    field.printToCSV(writeDirectory+"/displacementField.csv");
  }
  
  inline void DataRecord::writeColumnDisplacement(int bins) const {
    vector<pair<RealType, RealType> > columns[bins];
    Bounds sb = simData->getSimBounds();
    RealType binWidth = (sb.right-sb.left)/bins;
    // Bin data
    int domain_end = simData->getDomainEnd();
    RealType *px = simData->getPxPtr(), *py = simData->getPyPtr();
    int *it = simData->getItPtr();
    for (int i=0; i<domain_end; ++i) {
      if (-1<it[i]) {
	int b = (initialPositions.at(i).x-sb.left)/binWidth;
	auto point = pair<RealType, RealType>( initialPositions.at(i).y, py[i] );
	columns[b].push_back(point);
      }
    }

    // Write data
    string dirname = writeDirectory+"/Columns";
    mkdir(dirname.c_str(), 0777);
    for (int i=0; i<bins; ++i) printToCSV(dirname+"/col"+toStr(i)+".csv", columns[i]);
  }

  inline void DataRecord::getPressureData(const Bounds& region, ScalarField& field, RealType resolution) const {
    /*
    // Get the particle pressure data
    vector<PData> positions;
    sectors->getPressureData(positions); 
    // Find which particles cover which parts of the field
    int nsx = (region.right - region.left)/resolution, nsy = (region.top - region.bottom)/resolution;
    RealType sx = (region.right - region.left)/nsx, sy = (region.top - region.bottom)/nsy;
    ++nsx;
    ++nsy;
    vector<int> *sectors = new vector<int>[nsx*nsy];
    // Sectorize
    auto allParticles = simData->getParticles();
    for (const auto &p : allParticles) {
      int sec_x = (p.position.x - region.left)/sx;
      int sec_y = (p.position.y - region.bottom)/sy;
      if (-1<sec_x && sec_x<nsx && -1<sec_y && sec_y<nsy)
        sectors[nsx*sec_y+sec_x].push_back(i);
      // If part of the particle sticks into the region, place the particle in the closest sector. This ignores particles "diagonal" to the region that stick into the region
      else if (region.left<p.position.x+p.sigma   && -1<sec_y && sec_y<nsy)
        sectors[nsx*sec_y+0].push_back(i); // Left
      else if (p.position.x-p.sigma<region.right  && -1<sec_y && sec_y<nsy)
        sectors[nsx*sec_y+(nsx-1)].push_back(i); // Right
      else if (region.bottom<p.position.y+p.sigma && -1<sec_x && sec_x<nsx)
        sectors[sec_x].push_back(i); // Bottom
      else if (p.position.y-p.sigma<region.top    && -1<sec_x && sec_x<nsx)
        sectors[nsx*(nsy-1)+sec_x].push_back(i); // Top
      ++i;
    }

    // Set up field if that hasn't already been done
    if (field.empty()) {
      RealType sx = (region.right - region.left)/nsx, sy = (region.top - region.bottom)/nsy;
      field.setBounds(region);
      field.setResolution(sx);
      field.setPrintPoints(200);
    }
    
    RealType dx = field.getDx(), dy = field.getDy();
    // Update field
    for (int y=0; y<nsy; ++y)
      for (int x=0; x<nsx; ++x) {
	RealType p(0);
	int c(0);
	for (auto i : sectors[nsx*y+x]) {
	  p += std::get<5>(positions.at(i));
	  ++c;
	}
	if (c>0) field.at(x,y) += p/c;
      }
    
    for (const auto &p : positions) {
      int sec_x = (std::get<0>(p) - region.left)/sx;
      int sec_y = (std::get<1>(p) - region.bottom)/sy;
      if (-1<sec_x && sec_x<nsx && -1<sec_y && sec_y<nsy)
	field.at(sec_x, sec_y) += std::get<5>(p);
    }
    */
  }
  
  inline void DataRecord::getVortexData(const Bounds& region, ScalarField& field, RealType resolution, RealType cutoff) const {
    /*
    int nsx = (region.right - region.left)/resolution, nsy = (region.top - region.bottom)/resolution;
    RealType sdx = (region.right - region.left)/nsx, sdy = (region.top - region.bottom)/nsy;
    ++nsx;
    ++nsy;
    // Set up field if that hasn't already been done
    if (field.empty()) {
      RealType sx = (region.right - region.left)/nsx, sy = (region.top - region.bottom)/nsy;
      field.setBounds(region);
      field.setResolution(sx);
      field.setPrintPoints(200);
    }
    // Measure local angular momentum
    RealType X(0), Y(region.bottom + 0.5*sdy);
    for (int y=0; y<nsy; ++y, Y += sdy) {
      X = region.left + 0.5*sdx;
      for (int x=0; x<nsx; ++x, X += sdx) {
	vec2 pos(X, Y);
	auto particles = simData->getParticles(pos, cutoff);
	RealType ang = 0;
	for (const auto p : particles) {
	  if (p.sigma < 0.1) { // Radius restiction
	    vec2 displacement = simData->getDisplacement(p.position, pos);
	    RealType angM = displacement^((1./p.invMass)*p.velocity);
	    ang += angM;
	  }
	}
	field.at(x,y) = ang;
      }
    }
    */
  }

}
