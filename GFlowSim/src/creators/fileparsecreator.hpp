#ifndef __FILE_PARSE_CREATOR_HPP__GFLOW__
#define __FILE_PARSE_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"
#include "../utility/parsehelper.hpp"
#include "particletemplate.hpp"

#include <map>

namespace GFlowSimulation {

  const string Dimensions_Token = "Dimensions";
  const string Bounds_Token = "Bounds";
  const string Integrator_Token = "Integrator";
  const string Types_Token = "NTypes";
  const string Interactions_Token = "Force-grid";
  const string Boundary_Token = "Boundary-conditions";
  const string Fill_Token = "Fill-area";

  struct FillBounds {
    FillBounds(int dim) : bnd_dimensions(dim) {};

    //! @brief Returns the volume of the bounds
    virtual double vol()=0;
    //! @brief Return a random position within the bounds
    virtual void pick_position(RealType*)=0;
    //! @brief Get rectangular bounds that enclose the bounds
    virtual Bounds getBounds()=0;

    //! @brief Number of dimensions the bounds takes up.
    //!
    //! We include this so we can make lower dimensional fill areas - walls, lines, etc
    const int bnd_dimensions;
  };

  struct RectangularBounds : public FillBounds {
    RectangularBounds(int dim) : FillBounds(dim) {
      min = new RealType[dim];
      max = new RealType[dim];
    }

    RectangularBounds(Bounds bnds, int dim) : FillBounds(dim) {
      min = new RealType[dim];
      max = new RealType[dim];
      for (int d=0; d<dim; ++d) {
        min[d] = bnds.min[d];
        max[d] = bnds.max[d];
      }
    }

    ~RectangularBounds() {
      if (min) delete [] min;
      if (max) delete [] max;
    }

    virtual double vol() override {
      float v = 1;
      for (int d=0; d<bnd_dimensions; ++d)
        v *= (max[d] - min[d]);
      return v;
    }

    virtual void pick_position(RealType *x) override {
      for (int d=0; d<bnd_dimensions; ++d) 
        x[d] = drand48()*(max[d] - min[d]) + min[d];
    }

    Bounds getBounds() override {
      Bounds bnds(bnd_dimensions);
      for (int d=0; d<bnd_dimensions; ++d) {
        bnds.min[d] = min[d];
        bnds.max[d] = max[d];
      }
      return bnds;
    }

    float *min, *max;
  };

  struct SphericalBounds : public FillBounds {
    SphericalBounds(int dim) : FillBounds(dim), radius(0) {
      center = new RealType[dim];
    }

    ~SphericalBounds() {
      delete [] center;
    }

    virtual double vol() override {
      return sphere_volume(radius, bnd_dimensions);
    }

    virtual void pick_position(RealType *x) override {
      if (radius==0) {
        zeroVec(x, bnd_dimensions);
        return;
      }
      // Do this the dumb way for now, so it works in arbitrary # of dimensions
      bool good = false;
      while (!good) {
        for (int d=0; d<bnd_dimensions; ++d)
          x[d] = 2*(drand48() - 0.5)*radius + center[d];
        // Check whether the point is good
        good = sqr(x, bnd_dimensions)<=sqr(radius);
      }
    }

    Bounds getBounds() override {
      Bounds bnds(bnd_dimensions);
      for (int d=0; d<bnd_dimensions; ++d) {
        bnds.min[d] = center[d] - radius;
        bnds.max[d] = center[d] + radius;
      }
      return bnds;
    }

    float *center;
    float radius;
  };

  /**
  *  @brief A creator that creates a simulation from a file.
  *
  *  This creator parses the contents of a file and creates a simulation from the
  *  options specified in that file.
  *
  *  @see Creator
  */
  class FileParseCreator : public Creator {
  public:
    //! Constructor.
    FileParseCreator(int, char**);

    //! Constructor -- pass in a pointer to an ArgParse object.
    FileParseCreator(ArgParse*);

    //! Constructor -- pass in a pointer to an ArgParse object, and the configuration file name.
    FileParseCreator(ArgParse*, string);

    //! Create a simulation.
    virtual GFlow* createSimulation();

  private:

    inline void createFromOptions(HeadNode*);

    // --- Object selection

    inline Integrator* choose_integrator(HeadNode*) const;

    inline Interaction* choose_interaction(HeadNode*) const;

    inline BCFlag choose_bc(string&) const;

    inline void add_modifier(HeadNode*) const;

    // --- Creation helpers

    inline void fillArea(HeadNode*) const;

    inline FillBounds* getFillBounds(HeadNode*) const;

    inline void createParticle(HeadNode*) const;

    //! @brief Get all the headers with a certain heading, put into the supplied vector.
    inline void getAllMatches(string, vector<HeadNode*>&, std::multimap<string, HeadNode*>&) const;

    inline void getParticleTemplate(HeadNode*, std::map<string, ParticleTemplate>&) const;

    inline void makeRandomForces();

    inline RandomEngine* getRandomEngine(HeadNode*, string&) const;

    inline string copyFile() const;

    //! @brief The name of the file to load from
    string configFile;

    //! @brief The GFlow object the creator is creating.
    GFlow *gflow;

    //! @brief The number of particle types in the simulation
    int NTypes;

    //! @brief Particle templates usable anywhere in the configuration file.
    std::map<string, ParticleTemplate> global_templates;

    //! @brief The message the parser writes as it parses the configuration file
    mutable string parse_message, build_message;

    // Normal distribution
    mutable std::mt19937 generator;
    mutable std::normal_distribution<RealType> normal_dist;
  };

}
#endif // __FILE_PARSE_CREATOR_HPP__GFLOW__
