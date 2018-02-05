#include "FractureSimulation.hpp"
#include "../forces/TemperatureForce.hpp"
#include "../forces/ViscousDrag.hpp"

namespace GFlow {

  FractureSimulation::FractureSimulation() : SimulationBase(), limitTime(10), delay(0.002), bondRadius(0.5), markBreaks(true), lastBreak(0), breakDisplayLength(0.25), strengthen(true), breakDataSlot(-1), heatIters(0), maxIters(10), heatTime(5.), relaxTime(5.), heatTemperature(1.), heatRadius(0.5) {};

  void FractureSimulation::setUp(int argc, char** argv) {
    // Set the command line input
    SimulationBase::setUp(argc, argv);
    // Create data record, SimData, integrator
    dataRecord = new DataRecord;
    dataRecord->setCommand(argc, argv);
    breakDataSlot = dataRecord->getStatDataSlot("breaks times");
  }

  // Set parameters from the command line
  void FractureSimulation::parse() {
    // How long should the simulation run
    parser.get("time", limitTime);
    parser.get("markBreaks", markBreaks);
    parser.get("strengthen", strengthen);
    parser.get("maxIters", maxIters);
    parser.get("heatTime", heatTime);
    parser.get("relaxTime", relaxTime);
    parser.get("checkDelay", delay);
    parser.get("heatTemperature", heatTemperature);
    parser.get("heatRadius", heatRadius);

    // Standard parsing options and checking
    standardParsing();
  }

  // Run the simulation
  void FractureSimulation::run() {
    if (integrator==nullptr) return;

    // Set timers/counters
    bool running = true;
    RealType time = 0, timer = 0;
    heatIters = 0;

    // Initialize
    integrator->initialize(limitTime);

    // Set up verlet list tracking (do after integrator initialization)
    int num = simData->getDomainSize();
    bondLists.initialize(simData);
    
    // Find the average radius of particles
    RealType aveSg = 0.;
    for (int i=0; i<num; ++i) aveSg += simData->getSgPtr() [i];
    aveSg /= num;

    // Set radii and array sized
    bondRadius = 0.4*aveSg;
    bondLists.setSkinDepth(0);    
    bondLists.createVerletLists(); // Create the initial verlet lists
    verletListRecord.resize(num);  // Allocate enough space for each particle
    bondedParticles.resize(num);   // Allocate enough space for each particle
    accessor.resize(num);
    setVerletListRecord();
    
    // Make particles stationary
    simData->zeroMotion();

    //**
    // Add a moving wall
    Wall w(-2, 24, 2, 24);
    w.im = 1;
    w.vy = -1;
    w.rp = 0.288884;
    simData->addWall(w);
    // Add a bottom wall
    Wall w2(-2, 0, 2, 0);
    w2.rp = 0.288884;
    simData->addWall(w);

    // Pre
    integrator->preIntegrate();
    // Run the simulation
    while (running) {
      // Perform an integration step
      integrator->step();

      // Update the time and running flag
      time = integrator->getTime();
      running = integrator->isRunning();

      // Check for fragmentation - just not every timestep
      if (time - timer > delay) {
	// Check for breaks
	checkForBreaks();
	// Reset timer
	timer = time;
      }

      // If there is a break, reset and heat up the break
      if (strengthen && !breakages.empty()) {
	// This is the break time for the [heatIters] iteration
	if (-1<breakDataSlot) dataRecord->addStatData(breakDataSlot, RPair(heatIters, time));
	// Do a finite number of strengthening iterations
	if (heatIters < maxIters) {
	  // Print a message
	  if (!quiet) cout << "Iteration " << heatIters << ": t=" << simData->getTime() << ", " << numBreaks << " breaks.\n";
	  // Do the strengthening procedure
	  strengthenProcedure();
	  // Print a message
	  if (!quiet) cout << "Starting next trial.\n";
	  // Remake bond lists
	  setVerletListRecord();
	  // Reset the integrator time limit
	  integrator->setRunTime(limitTime);
	  // Reset timer/time
	  timer = 0;
	  time = 0;
	  // Reset break/bond count
	  numBreaks = 0;
	  numBonds = 0;
	}
	// We have done the requisite number of iterations
	else running = false;
      }
    }

    // Post integration
    integrator->postIntegrate();
  }

  // Save data from the simulation
  void FractureSimulation::write() {
    // Post - only print a ratio/runtime message if not doing strengthening
    standardWriting(!strengthen);
  }

  inline void FractureSimulation::checkForBreaks() {
    // Checks
    if (integrator==nullptr || integrator->getSectorization()==nullptr) return;
    
    // Update the verlet lists and get them
    bondLists.createVerletLists();
    auto& verletList = bondLists.getVerletList();

    // Clear the marked set
    marked.clear();

    // Set up the accesser
    for (auto &a : accessor) a = nullptr;
    for (auto &vl : verletList) accessor.at(*vl.begin()) = &vl;

    // Keep track of new bondings and breaking bonds
    vector<pair<int, int> > bonds, breaks;

    // Look for bonds and breaks
    compare(bonds, breaks);
    
    // Set verlet list record to be the current verlet list. We have to sort the data
    setVerletListRecord(verletList);

    // Get the current time, for recording when this happened
    RealType currentTime = integrator->getTime();
    // If there is bond data to record
    if (!bonds.empty()) {
      // Add bond pairs
      for (auto &b : bonds) {
	bondings.push_back(std::tuple<RealType,int,int>(currentTime, b.first, b.second));
	// We keep an unordered set of all the particles bonding or undergoing breakages
	if (!markBreaks) {
	  marked.insert(b.first);
	  marked.insert(b.second);
	}
      }
      // Update cumulative counting statistics
      numBonds += bonds.size();
      cumNumBonds.push_back(pair<RealType,int>(currentTime, numBonds));
    }
    // If there is breakage data to record
    if (!breaks.empty()) {
      // Add break pairs
      for (auto &b : breaks) {
	breakages.push_back(std::tuple<RealType,int,int>(currentTime, b.first, b.second));
	// We keep an unordered set of all the particles bonding or undergoing breakages
	if (markBreaks) {
	  marked.insert(b.first);
	  marked.insert(b.second);
	}
      }
      // Update cumulative counting statistics
      numBreaks += breaks.size();
      cumNumBreaks.push_back(pair<RealType,int>(currentTime, numBreaks));
      // Update last break time
      lastBreak = currentTime;
    }

    // Give the data record the marked particles
    if (!marked.empty() || currentTime-lastBreak>breakDisplayLength) dataRecord->setMarked(marked);
  }

  // Check for bonds and breaks
  inline void FractureSimulation::compare(vector<pair<int,int> >& bonds, vector<pair<int,int> >& breaks) {
    
    // If a formerly bonded particle went beyond a certain distance, it is no longer bonded
    for (int i=0; i<bondedParticles.size(); ++i) {
      for (auto j : bondedParticles.at(i))
	if (getHardCoreDistance(i,j) > 2*bondRadius)
	  breaks.push_back(pair<int,int>(i,j));
    }
    
    // Look for bonds
    // <--------------------
  }

  // Checks for bonds and breaks by comparing current and past verlet lists
  /*
  inline void FractureSimulation::compare(vector<pair<int,int> >& bonds, vector<pair<int,int> >& breaks) {
    // Lambda to find an element in a list
    auto check = [] (const int n, const VListSubType *list)->bool {
      if (list==nullptr) return false;
      for (auto &m : *list) if (m==n) return true;
      return false;
    };

    // Look for breaks
    for (int i=0; i<verletListRecord.size(); ++i) {
      if (!verletListRecord.at(i).empty()) {
	// For all of [i]'s neighbors, check that they stay neighbors
	for (auto nb : verletListRecord.at(i)) {
	  // Doesn't matter if you are in your verlet list. This is not a bonding/breakage
	  if (nb==i) continue;
	  // Check if [nb] is still in [i]'s neighbor list this time, or if [i] is in [nb]'s list this time
	  if (!check(nb, accessor.at(i)) && 
	      !check(i, accessor.at(nb)))
	    breaks.push_back(pair<int,int>(i,nb));
	}
      }
    }

    // Look for bonds
    for (int i=0; i<accessor.size(); ++i) {
      if (accessor.at(i)!=nullptr) {
	// For all of [i]'s neighbors, check that they were already neighbors. If not, a bond has formed
	for (auto nb : *accessor.at(i))
	  // Check if [nb] was in [i]'s neighbor list last time, or if [i] was in [nb]'s neighbor list last time
	  if (!check(nb, &verletListRecord.at(i)) && 
	      !check(i, &verletListRecord.at(nb)))
	    bonds.push_back(pair<int,int>(i,nb));
      }
    }
  }
  */

  inline void FractureSimulation::setVerletListRecord(VListType& verletList) {
    // Clear the old data
    for (int i=0; i<verletListRecord.size(); ++i) verletListRecord.at(i).clear();

    // Clear accessor
    for (int i=0; i<accessor.size(); ++i) accessor.at(i) = nullptr;

    // Insert new data
    for (const auto &nl : verletList) {
      int head = *nl.begin();
      verletListRecord.at(head) = nl; // If this causes an error, then there are new particles
    }

    // Lambda for finding objects in a list
    auto contains = [&] (int n, int p) ->bool {
      return std::find(bondedParticles.at(p).begin(), 
		       bondedParticles.at(p).end(), 
		       n)!=bondedParticles.at(p).end();
    };

    // Update [bondedParticles]
    vector<int> temp;
    RealType max = 0, min = 1;
    for (int i=0; i<verletListRecord.size(); ++i) {
      // If bonded particle is still in the (larger radius) neighbor list
      for(auto n : verletListRecord.at(i)) {
	if (n==i) continue;
	// Particle n is in the (big) neighborhood of i. Check if it was in the bond neighborhood, or if it should be (if it's close enough)
	// The rmin points (hard centers) of the particle should be within [bondRadius] of each other for them to count as being bound
	RealType hcDist = getHardCoreDistance(n,i);
	if ((hcDist < 2*bondRadius && contains(n,i)) || hcDist < bondRadius)
	  temp.push_back(n);
      }
      // Set the new bonded particles list
      bondedParticles.at(i) = temp;
      // Clear temp array
      temp.clear();
    }
  }

  inline void FractureSimulation::setVerletListRecord() {
    // Do checks
    if (integrator==nullptr) return;
    // Set up verlet list recording
    auto& verletList = bondLists.getVerletList();
    // Set the verlet list record
    setVerletListRecord(verletList);
  }

  inline void FractureSimulation::strengthenProcedure() {
    // Reset the scenario
    integrator->reset();
    simData->zeroMotion(); // A temporary fix since we don't store initial velocities, so resetting the simData doesn't reset velocities
    dataRecord->resetTimers(); // If we don't do this, data records 'last *' timers will be in the future, and no data will be collected
    // Get the particle numbers of the first breakage
    auto data = *breakages.begin();
    int p1 = std::get<1>(data), p2 = std::get<2>(data);
    // Reset breakages and bondings lists
    breakages.clear();
    bondings.clear();
    // Find the position of the particle that broke
    vec2 pos1(simData->getPxPtr() [p1], simData->getPyPtr() [p1]);
    vec2 pos2(simData->getPxPtr() [p2], simData->getPyPtr() [p2]);
    // Turn off particle's characteristics except for [Fixed]
    // If we don't fix the CV/CF/etc. particles, they can suddenly change position
    auto& characteristics = simData->getCharacteristics();
    vector<pair<int, Characteristic*> > characteristicRecord;
    for (auto& pc : characteristics) {
      // Keep track of what the original characteristic was
      characteristicRecord.push_back(pc);
      // Replace with a [Fixed] characteristic
      pc.second = new Fixed(pc.first, simData);
    }
    // Take away the simulation's forces. We will return them after the strengthening procedure
    auto externalForces = simData->transferExternalForces();
    // Apply a temperature near that position - pos, radius, temp
    ExternalForce *forces[3];
    forces[0] = new TemperatureForce(pos1, heatRadius, heatTemperature);
    forces[1] = new TemperatureForce(pos2, heatRadius, heatTemperature);
    forces[2] = new ViscousDrag;
    simData->addExternalForce(forces[0]);
    simData->addExternalForce(forces[1]);
    simData->addExternalForce(forces[2]);
    // Run for some amount of time
    integrator->integrate(heatTime);
    // Remove the external forces
    forces[0]->setFinished(true);
    forces[1]->setFinished(true);
    forces[2]->setFinished(true);
    simData->clearFinishedForces();
    // Add a viscous drag
    auto vd = new ViscousDrag(1.);
    simData->addExternalForce(vd);
    // Let the system relax
    integrator->integrate(relaxTime);
    // Remove external forces
    vd->setFinished(true);
    simData->clearFinishedForces(); 
    // Reset the time of the integrator and set the positions as the new initial positions
    simData->setInitialPositions();
    simData->zeroMotion(); // A temporary fix since we don't store initial velocities, so resetting the simData doesn't reset velocities
    integrator->reset();
    dataRecord->resetTimers(); // If we don't do this, data records 'last *' timers will be in the future, and no data will be collected
    // Return particles to their original characteristics
    for(auto pc: characteristicRecord) {
      Characteristic* &c = characteristics.at(pc.first);
      delete c;
      c = pc.second;
    }
    simData->resetCharacteristics();
    // Give back the original external force after resetting them
    for (auto f : externalForces) f->reset();
    simData->addExternalForce(externalForces);
    // Update count
    ++heatIters;
  }

  inline RealType FractureSimulation::getHardCoreDistance(int i, int j) {
    return sqrt(sqr(simData->getDisplacement(i,j)))-0.4*(simData->getSgPtr() [i] + simData->getSgPtr() [j]);
  }
  
}
