#ifndef __FILE_PARSE_CREATOR_HPP__GFLOW__
#define __FILE_PARSE_CREATOR_HPP__GFLOW__

#include "../base/creator.hpp"
#include "../utility/fileparse.hpp"
#include "particletemplate.hpp"

#include <map>

namespace GFlowSimulation {

  struct Molecule {
    ~Molecule() {
      for (auto &p : position) if (p) delete [] p;
    }
    //! @brief The position of the particles.
    vector<RealType*> position;
    //! @brief The mass of the particles.
    vector<RealType> mass;
    //! @brief The sigma of the particles.
    vector<RealType> sigma;
    //! @brief The type of the particle.
    vector<int> type;
    //! @brief Bond information. Comes in pairs.
    vector<int> bond;
    vector<RealType> relaxedLength;
    vector<RealType> bondStrength;
    //! @brief Angle information. Comes in triples.
    vector<int> angle;
    vector<RealType> relaxedAngle;
    vector<RealType> angleStrength;

    void set_com_coordinates() {
      RealType displacement[DIMENSIONS], weighted[DIMENSIONS];
      RealType totalMass(0);
      for (int i=0; i<position.size(); ++i) {
        copyVec(position[i], weighted);
        scalarMultVec(mass[i], weighted);
        plusEqVec(displacement, weighted);
        totalMass += mass[i];
      }
      scalarMultVec(1./totalMass, displacement);
      // Shift by center of mass
      for (int i=0; i<position.size(); ++i) {
        minusEqVec(position[i], displacement);
      }
    }
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

    //! @brief Exception class.
    struct UnexpectedOption {};
    struct BadStructure {
      BadStructure() : message("") {};
      BadStructure(string mess) : message(mess) {};
      string message;
    };

  private:

    inline void createFromOptions(GFlow*, std::multimap<string, HeadNode*>&);

    // --- Object selection

    inline Integrator* choose_integrator(string&) const;

    inline Interaction* choose_interaction(string&) const;

    inline BCFlag choose_bc(string&) const;

    // --- Creation helpers

    inline void fillArea(HeadNode*) const;

    //! @brief Get all the headers with a certain heading, put into the supplied vector.
    inline void getAllMatches(string, vector<HeadNode*>&, std::multimap<string, HeadNode*>&) const;

    inline void getParticleTemplate(HeadNode*, std::map<string, ParticleTemplate>&) const;

    inline RandomEngine* getRandomEngine(HeadNode*, string&) const;

    //! @brief The name of the file to load from
    string configFile;

    //! @brief The GFlow object the creator is creating.
    GFlow *gflow;

    //! @brief The bounds of the simulation we are creating
    Bounds bounds;

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
