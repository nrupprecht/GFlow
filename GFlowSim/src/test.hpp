#ifndef __TEST_HPP__GFLOW__
#define __TEST_HPP__GFLOW__

namespace GFlowSimulation {

  class Simulation {
    Simulation();

    void run(float);

  private:

    int number = 0;

    float *data = nullptr;

  };


}

#endif // __TEST_HPP__GFLOW__
