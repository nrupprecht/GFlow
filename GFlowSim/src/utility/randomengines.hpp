#ifndef __RANDOM_ENGINES_HPP__GFLOW__
#define __RANDOM_ENGINES_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  class RandomEngine {
  public:
    virtual ~RandomEngine() {};

    // \brief Generate a random number from the engine.
    virtual RealType generate() = 0;

    virtual RandomEngine* copy()=0;
  };

  class UniformRandomEngine : public RandomEngine {
  public:
    UniformRandomEngine(RealType l, RealType u) : lower(l), upper(u) {
      uniform_dist = std::uniform_real_distribution<RealType>(l, u);
    };
    //! \brief Generate a number.
    virtual RealType generate() override {
      return uniform_dist(global_generator);
    }

    virtual RandomEngine* copy() override {
      return new UniformRandomEngine(lower, upper);
    }
  private:
    RealType lower, upper;
    std::uniform_real_distribution<RealType> uniform_dist;
  };

  class ProportionalRandomEngine : public RandomEngine {
  public:
    ProportionalRandomEngine(RealType l, RealType u, int dims) : lower(l), upper(u), dimensions(dims) {
      uniform_dist = std::uniform_real_distribution<RealType>(0., 1.);
    };
    //! \brief Generate a number.
    virtual RealType generate() override {
      RealType x = uniform_dist(global_generator);

      // \todo This is only correct for 2 dimensions.
      
      return 1./(1./lower + (1./upper - 1./lower)*(1.-x));
    }

    virtual RandomEngine* copy() override {
      return new ProportionalRandomEngine(lower, upper, dimensions);
    }
  private:
    RealType lower, upper;
    int dimensions;
    std::uniform_real_distribution<RealType> uniform_dist;
  };

  class NormalRandomEngine : public RandomEngine {
  public:
    NormalRandomEngine(RealType a, RealType v) : ave(a), var(v) {
      nor_dist = std::normal_distribution<RealType>(a, v);
    };
    //! \brief Generate a number.
    virtual RealType generate() override {
      // Uses the global generator
      return nor_dist(global_generator);
    }
    virtual RandomEngine* copy() override {
      return new NormalRandomEngine(ave, var);
    }
  private:
    RealType ave, var;
    std::normal_distribution<RealType> nor_dist;
  };

  class DiscreteRandomEngine : public RandomEngine {
  public:
    DiscreteRandomEngine(vector<RealType>& prbs, vector<RealType>& vals) : probabilities(prbs), values(vals) {
      if (prbs.size()>vals.size()) throw false;
      discrete_dist = std::discrete_distribution<int>(prbs.begin(), prbs.end());
    };
    //! \brief Generate a number.
    virtual RealType generate() override {
      return values[discrete_dist(global_generator)];
    }
    virtual RandomEngine* copy() override {
      return new DiscreteRandomEngine(probabilities, values);
    }
  private:
    vector<RealType> probabilities, values;
    std::discrete_distribution<int> discrete_dist;
  };

  class DeterministicEngine : public RandomEngine {
  public:
    DeterministicEngine(RealType v) : value(v) {};
    //! \brief Generate a number.
    virtual RealType generate() override { return value; }

    virtual RandomEngine* copy() override {
      return new DeterministicEngine(value);
    }
  private:
    RealType value;
  };

}
#endif // __RANDOM_ENGINES_HPP__GFLOW__
