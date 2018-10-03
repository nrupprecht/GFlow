#ifndef __LINE_CRYSTAL_CREATOR_HPP__GFLOW__
#define __LINE_CRYSTAL_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"

namespace GFlowSimulation {

  //! @todo Commentary.
  class LineCrystalCreator : public Creator {
  public:
    //! Constructor.
    LineCrystalCreator(int, char**);

    //! Constructor -- pass in a pointer to an ArgParse object.
    LineCrystalCreator(ArgParse*);

    //! Seed random number generators.
    virtual void seedGenerator(uint);

    //! Create a simulation.
    virtual GFlow* createSimulation();

  private:

    //! @brief The structure function
    RealType gamma(RealType);

    //! @brief Radius of a sphere.
    RealType radius;

    //! @brief Distance between lines
    RealType lambda;

    //! @brief The width of the simulation - a parameter. 
    RealType width;

    //! @brief The height of the simulation - calculated as M*lambda.
    RealType height;

    //! @brief Number of spheres per line.
    int N;
    //! @brief Number of lines.
    int M;

    // Normal distribution
    std::mt19937 generator;
    std::normal_distribution<RealType> normal_dist;
  };

}
#endif // __LINE_CRYSTAL_CREATOR_HPP__GFLOW__