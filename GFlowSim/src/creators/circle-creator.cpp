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

    // Group net force object
    //shared_ptr<GroupNetForce> netforce = nullptr;
    if (track) {
      netforce = make_shared<GroupNetForce>(gflow);
      gflow->addDataObject(netforce);
    }

    // Place particles
    auto simData = gflow->getSimData(); // Get the simdata
    Vec pos(sim_dimensions), Zero(sim_dimensions);

    recursive_circle(0, sim_dimensions, center, pos, radius, process_bounds, simData);

    // if (sim_dimensions==2) {
    //   for (int i=0; i<n_particles; ++i, theta += dtheta) {
    //     // Position of next wall particles
    //     pos[0] = center[0] + radius*sin(theta);
    //     pos[1] = center[1] + radius*cos(theta);

    //     // If particle is in bounds, place it.
    //     if (process_bounds.contains(pos)) {
    //       // Get next global id
    //       int gid = simData->getNextGlobalID();
    //       // Add the particle to simdata. It is a particle of infinite mass.
    //       simData->addParticle(pos.data, Zero.data, sigma, 0, type);
    //       // Add particle to group net force data object.
    //       if (track) netforce->add(gid);
    //     }
    //   }
    // }
    // else if (sim_dimensions==3) {
    //   for (int i=0; i<n_particles && theta <= PI; ++i, theta += dtheta) {
    //     // Position of next ring of wall particles
    //     pos[0] = center[0] + radius*cos(theta);

    //     // Number of particles in the ring.
    //     real radius2 = radius*sin(theta);
    //     int n_particles_2 = max(static_cast<int>(ceil(0.5 * scale * (2*PI*radius2) / sigma)), 1);
    //     real dtheta2 = 2.*PI / n_particles_2, theta_2 = 0;

    //     for (int j=0; j<n_particles_2; ++j, theta_2 += dtheta2) {
    //       pos[1] = center[1] + radius2*cos(theta_2);
    //       pos[2] = center[2] + radius2*sin(theta_2);

    //       // If particle is in bounds, place it.
    //       if (process_bounds.contains(pos)) {
    //         // Get next global id
    //         int gid = simData->getNextGlobalID();
    //         // Add the particle to simdata. It is a particle of infinite mass.
    //         simData->addParticle(pos.data, Zero.data, sigma, 0, type);
    //         // Add particle to group net force data object.
    //         if (track) netforce->add(gid);
    //       }
    //     }
    //   }
    // }

  }

  void CircleCreator::recursive_circle(const int d, const int df, const Vec& center, Vec& pos, const real R, const Bounds& process_bounds, shared_ptr<SimData> simData) const {
    int n_particles = max(static_cast<int>(ceil(0.5 * scale * 2*PI*R / sigma)), 1);
    real dtheta = 2.*PI / n_particles;

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
          simData->addParticle(pos.data, nullptr, sigma, 0, type);
          // Add particle to group net force data object.
          if (track) netforce->add(gid);
        }
      }
    }
  }

}
