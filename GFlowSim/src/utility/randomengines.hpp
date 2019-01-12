#ifndef __RANDOM_ENGINES_HPP__GFLOW__
#define __RANDOM_ENGINES_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  class RandomEngine {
  public:
    virtual ~RandomEngine() {};

    // @brief Generate a random number from the engine.
    virtual double generate() = 0;

    virtual RandomEngine* copy()=0;
  };

  class UniformRandomEngine : public RandomEngine {
  public:
    UniformRandomEngine(double l, double u) : lower(l), upper(u) {
      uniform_dist = std::uniform_real_distribution<double>(l, u);
    };
    //! @brief Generate a number.
    virtual double generate() {
      return uniform_dist(global_generator);
    }

    virtual RandomEngine* copy() {
      return new UniformRandomEngine(lower, upper);
    }
  private:
    double lower, upper;
    std::uniform_real_distribution<double> uniform_dist;
  };

  class NormalRandomEngine : public RandomEngine {
  public:
    NormalRandomEngine(double a, double v) : ave(a), var(v) {
      nor_dist = std::normal_distribution<double>(a, v);
    };
    //! @brief Generate a number.
    virtual double generate() {
      // Uses the global generator
      return nor_dist(global_generator);
    }
    virtual RandomEngine* copy() {
      return new NormalRandomEngine(ave, var);
    }
  private:
    double ave, var;
    std::normal_distribution<double> nor_dist;
  };

  class DiscreteRandomEngine : public RandomEngine {
  public:
    DiscreteRandomEngine(vector<double>& prbs, vector<double>& vals) : probabilities(prbs), values(vals) {
      if (prbs.size()>vals.size()) throw false;
      discrete_dist = std::discrete_distribution<int>(prbs.begin(), prbs.end());
    };
    //! @brief Generate a number.
    virtual double generate() {
      return values[discrete_dist(global_generator)];
    }
    virtual RandomEngine* copy() {
      return new DiscreteRandomEngine(probabilities, values);
    }
  private:
    vector<double> probabilities, values;
    std::discrete_distribution<int> discrete_dist;
  };

  class DeterministicEngine : public RandomEngine {
  public:
    DeterministicEngine(double v) : value(v) {};
    //! @brief Generate a number.
    virtual double generate() { return value; }

    virtual RandomEngine* copy() {
      return new DeterministicEngine(value);
    }
  private:
    double value;
  };

}
#endif // __RANDOM_ENGINES_HPP__GFLOW__
