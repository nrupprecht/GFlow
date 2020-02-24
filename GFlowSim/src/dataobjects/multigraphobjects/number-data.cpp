#include "number-data.hpp"

namespace GFlowSimulation {

  NumberData::NumberData(GFlow *gflow) : MultiGraphObject(gflow, "NumberData", "time", "counts", gflow->getSimData()->ntypes()) {}; 

  void NumberData::post_step() {
    // Only record if enough time has gone by
    if (!DataObject::_check()) return;

    int ntypes = simData->ntypes();
    Vec count_population(ntypes);

    int size = simData->size_owned();
    for (int n=0; n<size; ++n) {
      int type = simData->Type(n);
      if (-1<type && type<ntypes)
        ++count_population[type];
    }

    // Create an entry with the average data. This function handles multiprocessor runs correctly.
    gatherData(gflow->getElapsedTime(), count_population);
  }
  

}