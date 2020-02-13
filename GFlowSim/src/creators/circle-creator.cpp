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
    // Check headings for validity.
    parser.check();

    // --- Parameters
    RealType radius = 1., sigma = 0.05;
    int type = 0;
    bool track = true;
    // --- Gather parameters
    Vec center = parser.argVec("Center");
    // Check that vectors are good.
    if (center.size()!=sim_dimensions) throw BadDimension("Center vector for Fill: Circle needs the correct dimensionality");
    // Get other data.
    parser.firstArg("Radius", radius);
    parser.firstArg("Sigma", sigma);
    parser.firstArg("Type", type);
    parser.firstArg("Track", track);

    // Compute number of particles, angle.
    RealType circumference = 2*PI*radius;
    int n_particles = max(static_cast<int>(ceil(0.5 * circumference / sigma)), 1);
    RealType dtheta = 2.*PI / n_particles, theta = 0;
    auto process_bounds = gflow->getTopology()->getProcessBounds();

    // Group net force object
    shared_ptr<GroupNetForce> netforce = nullptr;
    if (track) {
      netforce = make_shared<GroupNetForce>(gflow);
      gflow->addDataObject(netforce);
    }

    // Place particles
    auto simData = gflow->getSimData(); // Get the simdata
    Vec pos(sim_dimensions), Zero(sim_dimensions);
    for (int i=0; i<n_particles; ++i, theta += dtheta) {
      // Position of next wall particles
      pos[0] = center[0] + radius*sin(theta);
      pos[1] = center[1] + radius*cos(theta);

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
