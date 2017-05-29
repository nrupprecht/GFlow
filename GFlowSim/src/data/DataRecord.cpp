#include "DataRecord.hpp"
#include "../integrators/Integrator.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"
#include "../forces/ExternalForce.hpp"

namespace GFlow {
  DataRecord::DataRecord() : writeDirectory("RunData"), delay(1./15.), lastRecord(-delay), recIter(0), nsx(-1), nsy(-1), sdx(-1), sdy(-1), cutoff(-1), skinDepth(-1), recPos(false), lowerSizeLimit(0), recOption(1), recBulk(false), recBulkOutline(false), recDisplacementField(false), trackDisplacement(false), recPressField(false), recPerf(false), statRecPerf(-1), recMvRatio(false), statRecMvRatio(-1), recDt(false), statRecDt(-1), recDelay(false), statRecDelay(-1), center(false) {
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

  void DataRecord::record(SimData* simData, RealType time) {
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

    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it = simData->getItPtr();
    int domain_size = simData->domain_size;

    // Initialize displacement tracking, if neccessary
    if (recIter==0) {
      initialPositions = vector<vec2>(domain_size);
      for (int i=0; i<domain_size; ++i) initialPositions.at(i) = vec2(px[i], py[i]);
    }

    // Record position data
    if (recPos) {
      vector<PData> positions;
      if (recOption==1) // Record pressure
	simData->getPressureData(positions, lowerSizeLimit);
      else if (recOption==2) // Record by particle id
	recordByNumber(simData, positions);
      else if (recOption==3) // Record by number of verlet lists the particle is in
	recordByVerletList(simData, positions);
      else if (recOption==4) // Record by velocity
	recordByVelocity(simData, positions);
      else { // recOption==0 or default. Record nothing extra
	for (int i=0; i<domain_size; ++i)
	  if (it[i]>-1 && lowerSizeLimit<sg[i]) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], 0));
      }
      positionRecord.push_back(positions);
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
      getBulkData(simData, region, volumes, bulkField, lines);
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
      getPressureData(simData, region, pressField);
    }

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
      printToCSV(writeDirectory+"/Pos/number.csv", vector<int>(1,recIter));
      // Write position data
      for (int i=0; i<positionRecord.size(); ++i)
	if (!printToCSV(writeDirectory+"/Pos/pos", positionRecord.at(i), i))
	  std::cerr << "Failed to print to [" << writeDirectory << "/Pos/pos" << i << "].\n";

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
    if (simData && trackDisplacement) writeDisplacementData(simData);
    if (simData && recDisplacementField) writeDisplacementField(simData);

    // Write bulk data
    if (recBulk) {
      bulkField.printToCSV(writeDirectory+"/bulkField.csv");
      printToCSV(writeDirectory+"/bulkVolumes.csv", bulkVolumes);
    }

    // Write pressure data
    if (recPressField) pressField.printToCSV(writeDirectory+"/pressField.csv", 1./static_cast<RealType>(recIter));
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

      // Print checks
      if (0<simData->domain_size) writeParticleChecks(fout, simData);
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
    if (simData->sectors) {
      fout << "Sectorization summary (as of end of simulation):\n";
      int nsx = simData->sectors->getNSX(), nsy = simData->sectors->getNSY();
      
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

  vector<pair<RealType, RealType> > DataRecord::getStatFunctionData(int index) const {
    if (index<0 || statFunctionData.size()<=index) throw BadStatFunction(index);
    return statFunctionData.at(index);
  }

  string DataRecord::getStatFunctionName(int index) const {
    if (index<0 || statFunctionData.size()<=index) throw BadStatFunction(index);
    return statFunctionName.at(index);
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

  void DataRecord::writeParticleChecks(std::ofstream& fout, SimData* simData) const {
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

  void DataRecord::recordByNumber(SimData* simData, vector<PData>& positions) const {
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
      if (it[i]>-1 && lowerSizeLimit<sg[i]) positions.push_back(pos.at(i));
  }

  void DataRecord::recordByVelocity(SimData* simData, vector<PData>& positions) const {
    // Get the arrays
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    RealType *vx = simData->getVxPtr();
    RealType *vy = simData->getVyPtr();
    int *it = simData->getItPtr();
    int domain_size = simData->domain_size;

    // We will sort out particles with it<0 at the end
    for (int i=0; i<domain_size; ++i)
      if (it[i]>-1 && lowerSizeLimit<sg[i]) positions.push_back(PData(px[i], py[i], sg[i], th[i], it[i], sqrt(sqr(vx[i])+sqr(vy[i]))));
  }

  void DataRecord::getBulkData(SimData* simData, const Bounds& region, vector<RealType>& volumes, ScalarField& field, vector<pair<vec2,vec2> >& lines, RealType resolution, RealType dr, RealType lowerVCut, RealType upperVCut) const {
    int nsx = (region.right - region.left)/resolution, nsy = (region.top - region.bottom)/resolution;
    RealType sx = (region.right - region.left)/nsx, sy = (region.top - region.bottom)/nsy;
    ++nsx;
    ++nsy;
    list<int> *sectors = new list<int>[nsx*nsy];
    // Fill sectors
    int i=0;
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
	for (const auto index : sectors[nsx*y+x]) {
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
	      for (const auto index : sectors[nsx*Y+X]) {
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
      delete [] sectors;
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
    delete [] sectors;
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

  inline void DataRecord::writeDisplacementData(SimData* simData) const {
    // Get data
    int domain_size = simData->getDomainSize();
    RealType *px = simData->getPxPtr();
    RealType *py = simData->getPyPtr();
    RealType *sg = simData->getSgPtr();
    RealType *th = simData->getThPtr();
    int *it      = simData->getItPtr();
    // We may omit the largest particle
    RealType maxSg = StatFunc_MaxR_Sigma(simData);
    // Collect data
    vector<PData> displacementI;
    vector<PData> displacementF;
    for (int i=0; i<domain_size; ++i) {
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

  inline void DataRecord::writeDisplacementField(SimData* simData) const {
    // Get data
    int domain_size = simData->getDomainSize();
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
    for (int i=0; i<domain_size; ++i) {
      RealType disp = sqrt(sqr(initialPositions.at(i) - vec2(px[i], py[i])));
      // Bin data
      int sec_x = (initialPositions.at(i).x-bounds.left)*idx;
      int sec_y = (initialPositions.at(i).y-bounds.bottom)*idy;
      ++count[sec_y*nsx+sec_x];
      field.at(sec_x, sec_y) += disp;
    }
    for (int y=0; y<nsy; ++y)
      for (int x=0; x<nsx; ++x)
	if (count[y*nsx+x]>0) field.at(x,y) /= (RealType)count[y*nsx+x];
    // Write displacement to file
    field.printToCSV(writeDirectory+"/displacementField.csv");
  }

  inline void DataRecord::getPressureData(SimData* simData, const Bounds& region, ScalarField& field, RealType resolution) const {
    // Get the particle pressure data
    vector<PData> positions;
    simData->getPressureData(positions); 
    // Find which particles cover which parts of the field
    int nsx = (region.right - region.left)/resolution, nsy = (region.top - region.bottom)/resolution;
    RealType sx = (region.right - region.left)/nsx, sy = (region.top - region.bottom)/nsy;
    ++nsx;
    ++nsy;
    /*
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
    */

    // Set up field if that hasn't already been done
    if (field.empty()) {
      RealType sx = (region.right - region.left)/nsx, sy = (region.top - region.bottom)/nsy;
      field.setBounds(region);
      field.setResolution(sx);
      field.setPrintPoints(200);
    }
    
    /*
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
    */
    
    for (const auto &p : positions) {
      int sec_x = (std::get<0>(p) - region.left)/sx;
      int sec_y = (std::get<1>(p) - region.bottom)/sy;
      if (-1<sec_x && sec_x<nsx && -1<sec_y && sec_y<nsy)
	field.at(sec_x, sec_y) += std::get<5>(p);
    }
    
  }

}
