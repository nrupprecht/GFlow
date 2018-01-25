#include "Creator.hpp"

namespace GFlow {
  
  void Creator::createBuoyancy(SimData *& simData, Integrator *& integrator, RealType radius, RealType density, vec2 velocity, bool constant, bool insert) {
    Bounds sb = simData->getSimBounds();
    // Change the bounds for the addition of the new particle
    RealType height = StatFunc_HighestBall(simData); // Height of the tops of the balls
    Bounds nb(sb.left, sb.right, sb.bottom, height+6*radius);
    // Set new bounds
    simData->setSimBounds(nb);
    simData->setBounds(nb);
    // Get rid of the old walls
    simData->getWalls().clear();
    // Create new walls
    simData->addWall(nb);
    // Insert the intruding particle
    if (radius>0) {
      Particle P(0.5*(nb.right+nb.left), height + radius, radius);
      P.setDensity(density);
      P.velocity = velocity;
      // Add either with inserter, constant velocity, or as normal
      if (insert) simData->addParticle(P, new Insertion(velocity, 0));
      else if (constant) simData->addParticle(P, new ConstantVelocity(velocity, 0, true));
      else simData->addParticle(P);
    }

    // Make sure we use a small time step
    if (integrator) reinterpret_cast<VelocityVerletIntegrator*>(integrator)->setMaxTimeStep(1e-4);
  }

  void Creator::createAero(SimData *& simData, Integrator *& integrator, RealType radius, RealType density, vec2 velocity, bool constant) {
    Bounds sb = simData->getSimBounds();
    // Change the bounds for the addition of the new particle
    RealType height = StatFunc_HighestBall(simData); // Height of the tops of the balls
    Bounds nb(sb.left, sb.right, sb.bottom, height);
    // Set new bounds
    simData->setSimBounds(nb);
    simData->setBounds(nb);
    // Get rid of the old walls
    simData->getWalls().clear();
    // Set wrapping in the y-direction
    simData->setWrapY(true);
    // Remove external forces
    simData->clearExternalForces();
    // Create new walls
    Wall lw(nb.left, nb.bottom, nb.left, nb.top);
    Wall rw(nb.right, nb.bottom, nb.right, nb.top);
    simData->addWall(lw); simData->addWall(rw);
    // Insert the intruding particle
    if (radius>0) {
      Particle P(0.5*(nb.right+nb.left), 0.5*(nb.top+nb.bottom), radius);
      P.setDensity(density);
      P.velocity = velocity;
      // Add either with constant velocity, or as normal
      if (constant) simData->addParticle(P, new ConstantVelocity(velocity, 0.));
      else simData->addParticle(P);
      // Make sure there is no serious overlap
      StandardSectorization remover(simData);
      remover.removeOverlapping();
    }
  }

  void Creator::createMixer(SimData *& simData, Integrator *& integrator, RealType radius, RealType circleRadius, RealType omega, bool gravity) {
    Bounds sb = simData->getSimBounds();
    // Change the bounds for the addition of the new particle
    RealType height = StatFunc_HighestBall(simData); // Height of the tops of the balls
    Bounds nb(sb.left, sb.right, sb.bottom, height);
    // Set new bounds
    simData->setSimBounds(nb);
    simData->setBounds(nb);
    // Get rid of the old walls
    simData->getWalls().clear();
    // Set wrapping in the x-direction
    simData->setWrapX(true);
    // Create new walls
    Wall bw(nb.left, nb.bottom, nb.right, nb.bottom);
    Wall tw(nb.left, nb.top, nb.right, nb.top);
    Wall lw(nb.left, nb.bottom, nb.left, nb.top);
    Wall rw(nb.right, nb.bottom, nb.right, nb.top);
    simData->addWall(bw); 
    simData->addWall(tw);
    simData->addWall(lw);
    simData->addWall(rw);
    if (!gravity) simData->clearExternalForces();
    // Insert the mixing particle
    Particle P(0.5*(nb.right+nb.left)+circleRadius, 0.5*(nb.top+nb.bottom), radius);
    P.setDensity(1.); // Density doesn't matter, its motion is controlled
    // Add either with constant velocity, or as normal
    simData->addParticle(P, new Circulate(circleRadius, omega));
    // Make sure there is no serious overlap
    StandardSectorization remover(simData);
    remover.removeOverlapping();
    // Make sure we use a small time step
    if (integrator) reinterpret_cast<VelocityVerletIntegrator*>(integrator)->setMaxTimeStep(1e-4);
  }

  bool Creator::createRegion(Region& region, SimData* simData) {
    // Make sure the region is within the bounds
    if (!simData->simBounds.contains(region.bounds)) return false;

    // List of particles we are creating
    vector<Particle> particles;

    // First create radii list (this also allocates the correct number of particles)
    if (region.sigma) region.sigma->makeValues(region, particles);
    else return false; // We require a sigma function

    // Assign interaction
    if (region.interaction) region.interaction->makeValues(region, particles);
    else Homogeneous_Interaction().makeValues(region, particles);

    // Assign positions
    if (region.position) region.position->makeValues(region, particles);
    else Uniform_Space_Distribution().makeValues(region, particles);
    
    // Assign inertias (masses and moments of inertia)
    if (region.inertia) region.inertia->makeValues(region, particles);
    else Constant_Density().makeValues(region, particles);

    // Assign repulsion
    if (region.repulsion) region.repulsion->makeValues(region, particles);
    else Uniform_Random_Repulsion().makeValues(region, particles);

    // Assign dissipation
    if (region.dissipation) region.dissipation->makeValues(region, particles);
    else Uniform_Random_Dissipation().makeValues(region, particles);

    // Assign coefficient of friction
    if (region.coeff) region.coeff->makeValues(region, particles);
    else Uniform_Random_Coeff().makeValues(region, particles);

    // Let the simulation relax, using heavy viscous drag, so particles do not overlap with things
    SimData *relax = new SimData(region.bounds, region.bounds);
    relax->setWrap(false);
    relax->reserve(particles.size());
    relax->addParticle(particles);
    // Add containment walls
    relax->addWall(region.bounds);
    // Add all walls from the simulation
    relax->addWall(simData->getWalls());
    // Create the relaxation integrator
    VelocityVerletIntegrator verlet(relax);
    verlet.setAdjustTimeStep(false);
    verlet.addExternalForce(new ViscousDrag(1.)); // This is large enough
    RealType phi = relax->getPhi();
    verlet.integrate(2*sqr(phi));
    // Remove drag force
    particles = relax->getParticles();
    delete relax;
    
    // Set velocities and omegas
    // Assign velocity
    if (region.velocity) region.velocity->makeValues(region, particles);
    else Normal_Random_Velocity().makeValues(region, particles);

    // Assign theta
    if (region.theta) region.theta->makeValues(region, particles);
    else Uniform_Random_Theta().makeValues(region, particles);

    // Assign omega
    if (region.omega) region.omega->makeValues(region, particles);
    else Normal_Random_Omega().makeValues(region, particles);    
    
    // Reserve more room for the new particles
    simData->reserveAdditional(particles.size(), 0);
    // Add the particles to simData
    simData->addParticle(particles);

    // Return success
    return true;
  }
  
}
