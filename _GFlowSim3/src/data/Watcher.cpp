#include "Watcher.hpp"
#include "../control/SimData.hpp"
#include "DataRecord.hpp"

#include <string>
using std::string;

namespace GFlow {

  Watcher::Watcher(int i, bool t, DataRecord *d, string tag) : id(i), type(t), dataRecord(d), slot(-1), name("") {
    // Create the name
    if (t) name += "particle_";
    else name += "wall_";
    name += toStr(i);
    name += ("_" + tag);
    // Get a slot
    if (dataRecord) slot = dataRecord->getStatDataSlot(name);
  }

  void Watcher::setDataRecord(DataRecord *d) {
    dataRecord = d;
    slot = dataRecord->getStatDataSlot(name);
  }

  SimData* Watcher::getSimData() {
    return dataRecord->simData;
  }

  NetForceWatcher::NetForceWatcher(int i, bool t, DataRecord *d) : Watcher(i, t, d, "net_force") {};

  void NetForceWatcher::record() {
    // Retrieve the SimData
    SimData *sd = getSimData();
    RealType f = 0;
    // Particle
    if (type) {
      vec2 F (sd->getFxPtr() [id], sd->getFyPtr() [id]);
      f = sqrt(sqr(f));
    }
    // Wall
    else {
      vec2 F (sd->getWalls() [id].fx, sd->getWalls() [id].fy);
      f = sqrt(sqr(F));
    }
    // Record the data
    dataRecord->addStatData(slot, f);
  }

  LKEWatcher::LKEWatcher(int i, bool t, DataRecord *d) : Watcher(i, t, d, "linear_ke") {};

  void LKEWatcher::record() {
    // Retrieve the SimData
    SimData *sd = getSimData();
    RealType ke = 0; 
    // Particle
    if (type)
      ke = 0.5*(1./sd->getImPtr() [id]) * 
      (sqr(sd->getVxPtr() [id]) + 
       sqr(sd->getVyPtr() [id]));
    // Wall
    else 
      ke = 0.5*(1./sd->getWalls() [id].im) *
	(sqr(sd->getWalls() [id].vx) +
	 sqr(sd->getWalls() [id].vy));
    // Record the data
    dataRecord->addStatData(slot, ke);
  }

  PyWatcher::PyWatcher(int i, bool t, DataRecord *d) : Watcher(i, t, d, "py") {}

  void PyWatcher::record() {
    // Retrieve the SimData
    SimData *sd = getSimData();
    RealType py = 0;
    // Particle
    if (type) py = sd->getPyPtr() [id];
    else      py = sd->getWalls() [id].py;
    // Record the data
    dataRecord->addStatData(slot, py);
  }

}
