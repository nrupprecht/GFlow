#include "wallcreator.hpp"
// Other files
#include "../utility/treeparser.hpp"
#include "../allmodifiers.hpp"

namespace GFlowSimulation {

  void WallCreator::createArea(HeadNode *head, GFlow *gflow, const std::map<string, string>& variables, vector<ParticleFixer>& particle_fixers) {
    // Check if head is good
    if (head==nullptr) return;
    // For local particle templates
    std::map<string, ParticleTemplate> particle_templates;

    // Get the dimensions
    int sim_dimensions = gflow->getSimDimensions();

    // Create parser, with variables.
    TreeParser parser(head, variables);
    // Declare valid options.
    parser.addHeadingNecessary("Start", "We need starting position!");
    parser.addHeadingNecessary("End", "We need ending position!");
    parser.addHeadingNecessary("Radius", "We need particle radii!");
    parser.addHeadingOptional("Spacing");
    parser.addHeadingOptional("Type");
    parser.addHeadingOptional("Modifier");
    // Check headings for validity.
    parser.check();

    // --- Parameters
    RealType radius = 0.05;
    RealType spacing = 0.05;
    int type = 0;

    // --- Gather parameters
    Vec start = parser.argVec("Start");
    Vec end   = parser.argVec("End");
    Vec norm  = end - start;
    norm.normalize();
    // Check that vectors are good.
    if (start.size()!=sim_dimensions) throw BadDimension("Start point vector for Fill: Wall needs the correct dimensionality");
    if (end.size()!=sim_dimensions)   throw BadDimension("End point vector for Fill: Wall needs the correct dimensionality");

    // Get other data.
    parser.firstArg("Radius", radius);
    parser.firstArg("Type", type);

    // Default value for spacing is radius.
    spacing = 2*radius;
    parser.firstArg("Spacing", spacing);
    
    // Correct spacing to fit the distance between starting and ending points.
    RealType distance = distanceVec(start, end);
    int n_particles = max(static_cast<int>(ceil((distance - 2*radius) / spacing)), 1);
    spacing = distance / n_particles;

    // Place particles
    SimData *simData = gflow->getSimData(); // Get the simdata
    RealType lambda = radius;
    Group wall_group;
    Vec pos(sim_dimensions), Zero(sim_dimensions);
    for (int i=0; i<n_particles; ++i) {
      // Position of next wall particles
      pos = start + lambda * norm;
      // Get next global id
      int gid = simData->getNextGlobalID();
      // Add the particle to simdata
      simData->addParticle(pos.data, Zero.data, radius, 0, type);
      // Add the particle to the group
      wall_group.add(gid);
      // Increment parameter
      lambda += spacing;
    }

    // If there is a modifier, set that up.
    if (parser.begin("Modifier")) {
      do {
        if (parser.argName()=="Demon") {
          // Create a demon
          Demon *demon = new Demon(gflow);
          Vec center(sim_dimensions);
          gflow->getBounds().center(center.data);
          // Set partition position center.
          demon->setPartitionPosition(center[0]);
          // Set group.
          *dynamic_cast<Group*>(demon) = wall_group;
          RealType tau = 0;
          if (parser.firstArg("Tau", tau)) demon->setTau(tau);
          // Add the demon
          gflow->addModifier(demon);
        }
        else throw false; // BadStructure("Unrecognized modifier in wall creator, ["+parser.argName()+"].");
      } while (parser.next());
    }

  }

}