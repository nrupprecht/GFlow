#include "Simulator.h"

int main() {

  Simulator simulation;

  simulation.setSectorize(true);
  cout << "Sectorize: True" << endl;

  int total = 30;
  cout << "{";
  for (int i=0; i<total; i++) {
    simulation.createHopper(10*(i+1));
    simulation.run(10);
    cout << "{" << (i+1)*10 << "," << simulation.getRunTime() << "}";
    if (i!=total-1) cout << ",";
  }
  cout << "}\n";
  
  return 0;
}

