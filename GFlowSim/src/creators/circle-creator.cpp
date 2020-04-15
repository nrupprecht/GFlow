#include "circle-creator.hpp"
// Other files


namespace GFlowSimulation {

  void CircleCreator::createArea(HeadNode *head, GFlow *gflow, const std::map<string, string>& variables, vector<ParticleFixer>& particle_fixers) {
    // Check if head is good
    if (head==nullptr) return;

    // Get the dimensions
    int sim_dimensions = gflow->getSimDimensions();

    // Create parser, with variables.
    TreeParser parser(head, variables);
    // Declare valid options.
    parser.addHeadingNecessary("Center", "We need circle center!");
    parser.addHeadingOptional("Radius"); // Circle radius
    parser.addHeadingOptional("Sigma");
    parser.addHeadingOptional("Type");
    parser.addHeadingOptional("Scale");
    // Check headings for validity.
    parser.check();

    // --- Gather parameters
    Vec center = parser.argVec("Center");
    // Check that vectors are good.
    if (center.size()!=sim_dimensions) throw BadDimension("Center vector for Fill: Circle needs the correct dimensionality");
    // Get other data.
    parser.firstArg("Radius", radius);
    parser.firstArg("Sigma", sigma);
    parser.firstArg("Type", type);
    parser.firstArg("Track", track);
    parser.firstArg("Scale", scale);

    // If the radius is zero or negative, that means that we should not actually create a circle. Return here.
    if (radius<=0) return;

    // Compute number of particles, angle.
    real circumference = 2*PI*radius;
    int n_particles = max(static_cast<int>(ceil(0.5 * scale * circumference / sigma)), 1);
    real dtheta = 2.*PI / n_particles, theta = 0;
    auto process_bounds = gflow->getTopology()->getProcessBounds();

    // Group net force object.
    if (track) {
      netforce = make_shared<GroupNetForce>(gflow);
      gflow->addDataObject(netforce);
    }

    auto simData = gflow->getSimData(); // Get the simdata.
    Vec pos(sim_dimensions), Zero(sim_dimensions);
    // Recursively create a sphere out of lower dimensional spheres.
    recursive_circle(0, sim_dimensions, center, pos, radius, process_bounds, simData);
  }

  void CircleCreator::recursive_circle(const int d, const int df, const Vec& center, Vec& pos, const real R, const Bounds& process_bounds, shared_ptr<SimData> simData) const {
    int n_particles = max(static_cast<int>(ceil(0.5 * scale * 2*PI*R / sigma)), 1);
    real dtheta = 2.*PI / n_particles;

    Vec Zero(df);
    bool final_circle = (d+2==df);
    for (real theta = 0.; theta<PI || (final_circle && theta<2*PI); theta += dtheta) {
      real Rprime = R*sin(theta);
      pos[d]   = center[d]   + R*cos(theta);
      pos[d+1] = center[d+1] + Rprime;
      // (Possible) recursive part.
      if (!final_circle) recursive_circle(d+1, df, center, pos, Rprime, process_bounds, simData);
      // Place a particle.
      else {
        // If particle is in bounds, place it.
        if (process_bounds.contains(pos)) {
          // Get next global id
          int gid = simData->getNextGlobalID();
          // Add the particle to simdata. It is a particle of infinite mass.
          simData->addParticle(pos.data, Zero.data, sigma, 0, type);
          // Add particle to group net force data object.
          if (track) netforce->add(gid);
        }
      }
    }
  }

}
