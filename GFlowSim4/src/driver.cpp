#include "gflow.hpp"

int main(int argc, char **argv) {

  GFlowSimulation::GFlow gflow(argc, argv);

  gflow.run(1.);

  return 0;
}
