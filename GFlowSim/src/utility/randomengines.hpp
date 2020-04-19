#ifndef __RANDOM_ENGINES_HPP__GFLOW__
#define __RANDOM_ENGINES_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  class RandomEngine {
  public:
    //! \brief Virtual destructor. Keeps warning messages from happening.
    virtual ~RandomEngine() {};

    // \brief Generate a random number from the engine.
    virtual real generate() const = 0;

    //! \brief Create a copy of this random engine.
    virtual RandomEngine* copy() const = 0;
  };

  class UniformRandomEngine : public RandomEngine {
  public:
    UniformRandomEngine(real l, real u) : lower(l), upper(u) {
      uniform_dist = std::uniform_real_distribution<real>(l, u);
    };
    //! \brief Generate a number.
    virtual real generate() const override {
      return uniform_dist(global_generator);
    }

    virtual RandomEngine* copy() const override {
      return new UniformRandomEngine(lower, upper);
    }
  private:
    real lower, upper;
    mutable std::uniform_real_distribution<real> uniform_dist;
  };

  class ProportionalRandomEngine : public RandomEngine {
  public:
    ProportionalRandomEngine(real l, real u, int dims) : lower(l), upper(u), dimensions(dims) {
      uniform_dist = std::uniform_real_distribution<real>(0., 1.);
    };
    //! \brief Generate a number.
    virtual real generate() const override {
      real x = uniform_dist(global_generator);
      // \todo This is only correct for 2 dimensions.
      return 1./(1./lower + (1./upper - 1./lower)*(1.-x));
    }

    virtual RandomEngine* copy() const override {
      return new ProportionalRandomEngine(lower, upper, dimensions);
    }
  private:
    real lower, upper;
    int dimensions;
    mutable std::uniform_real_distribution<real> uniform_dist;
  };

  class NormalRandomEngine : public RandomEngine {
  public:
    NormalRandomEngine(real a, real v) : ave(a), var(v) {
      nor_dist = std::normal_distribution<real>(a, v);
    };
    //! \brief Generate a number.
    virtual real generate() const override {
      // Uses the global generator
      return nor_dist(global_generator);
    }
    virtual RandomEngine* copy() const override {
      return new NormalRandomEngine(ave, var);
    }
  private:
    real ave, var;
    mutable std::normal_distribution<real> nor_dist;
  };

  class DiscreteRandomEngine : public RandomEngine {
  public:
    DiscreteRandomEngine(const vector<real>& prbs, const vector<real>& vals) : probabilities(prbs), values(vals) {
      if (prbs.size()>vals.size()) throw false;
      discrete_dist = std::discrete_distribution<int>(prbs.begin(), prbs.end());
    };
    //! \brief Generate a number.
    virtual real generate() const override {
      return values[discrete_dist(global_generator)];
    }
    virtual RandomEngine* copy() const override {
      return new DiscreteRandomEngine(probabilities, values);
    }
  private:
    vector<real> probabilities, values;
    mutable std::discrete_distribution<int> discrete_dist;
  };

  class DeterministicEngine : public RandomEngine {
  public:
    DeterministicEngine(real v) : value(v) {};
    //! \brief Generate a number.
    virtual real generate() const override { return value; }

    virtual RandomEngine* copy() const override {
      return new DeterministicEngine(value);
    }
  private:
    real value;
  };

}
#endif // __RANDOM_ENGINES_HPP__GFLOW__
