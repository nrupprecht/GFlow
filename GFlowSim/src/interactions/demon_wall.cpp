#include "demon_wall.hpp"
// Other files
#include "../base/interactionhandler.hpp"

namespace GFlowSimulation {

  DemonWall::DemonWall(GFlow *gflow) : Interaction(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION) {};

  void DemonWall::interact() const {
    // Do dimensional check.
    // \todo Should probably have some sort of global error message system.
    if (sim_dimensions!=2) return;

    // Get the data pointers.
    RealType **x = Base::simData->X();
    RealType **v = Base::simData->V();
    RealType **f = Base::simData->F();
    RealType *sg = Base::simData->Sg();
    RealType *im = Base::simData->Im();
    int    *type = Base::simData->Type();
    int    entry = simData->getIntegerData("Side");
    int    *side = nullptr;
    if (0<=entry) side = Base::simData->IntegerData(entry); // Left: 0, Right: 1

    // Make sure all needed pointers are non null.
    // \todo Should probably have some sort of global error message system.
    if (x==nullptr || v==nullptr || sg==nullptr || im==nullptr || type==nullptr || side==nullptr) return;

    // Get the bounds and boundary conditions
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[2];
    copyVec(Base::gflow->getBCs(), boundaryConditions, 2); // Keep a local copy of the bcs
    // Extract bounds related data
    RealType bnd_x = bounds.wd(0);
    RealType bnd_y = bounds.wd(1);

    // Needed constants
    RealType sg1, sg2, dx, dy, rsqr, r, invr, magnitude;
    Vec X(2);
    vector<int> neighbors;
    vector<int> moved;

    // --- Go through all particles
    for (int i=0; i<verlet.size(); i+=2) {
      int id1 = verlet[i];
      int id2 = verlet[i+1];
      // Check if the types are good
      if (type[id1]<0 || type[id2]<0) continue;
      // The object of larger type will be the wall, make this id2
      if (type[id1]>type[id2]) std::swap(id1, id2);
      // Calculate displacement.
      dx = x[id1][0] - x[id2][0];
      dy = x[id1][1] - x[id2][1];
      // Harmonic corrections to distance.
      if (boundaryConditions[0]==BCFlag::WRAP) {
        RealType dX = bnd_x - fabs(dx);
        if (dX<fabs(dx)) dx = dx>0 ? -dX : dX;
      }  
      if (boundaryConditions[1]==BCFlag::WRAP) {
        RealType dY = bnd_y - fabs(dy);
        if (dY<fabs(dy)) dy = dy>0 ? -dY : dY;
      }  
      // Calculate squared distance
      rsqr = dx*dx + dy*dy;
      // Get radii
      sg1 = sg[id1];
      sg2 = sg[id2];
      // If close, interact.
      if (0<rsqr && rsqr < sqr(sg1 + sg2)) {
        // If just turned on, remove all particles touching us.
        if (turn_on) {
          if (side[id1]==0 && im[id1]>0) {
            // Place the particle somewhere randomly on its side, where it does not overlap with anything else.
            do {
              neighbors.clear();
              X[0] = bounds.min[0] + 0.025*bnd_x + drand48()*0.45*bnd_x; // Make sure we don't appear in a wall.
              X[1] = bounds.min[1] + drand48()*bnd_y;
              handler->getAllWithin(X, neighbors);
            } while (!neighbors.empty());
            // Place particle.
            x[id1][0] = X[0];
            x[id1][1] = X[1];
            v[id1][0] = -v[id1][0];
            moved.push_back(id1);
          }
          else if (side[id1]==1 && im[id1]>0) {
            // Place the particle somewhere randomly on its side, where it does not overlap with anything else.
            do {
              neighbors.clear();
              X[0] = bounds.min[0] + 0.525*bnd_x + drand48()*0.45*bnd_x; // Make sure we don't appear in a wall.
              X[1] = bounds.min[1] + drand48()*bnd_y;
              handler->getAllWithin(X, neighbors);
            } while (!neighbors.empty());
            // Place particle.
            x[id1][0] = X[0];
            x[id1][1] = X[1];
            v[id1][0] = -v[id1][0];
            moved.push_back(id1);
          }
        }
        // Normal interaction
        else {
          // Calculate distance, inverse distance.
          r = sqrt(rsqr);
          invr = 1./r;
          // Create a normal vector
          dx *= invr;
          dy *= invr;
          // Calculate the magnitude of the force
          magnitude = repulsion*(sg1 + sg2 - r);
          // Update forces
          f[id1][0] += magnitude * dx;
          f[id2][0] -= magnitude * dx;
          f[id1][1] += magnitude * dy;
          f[id2][1] -= magnitude * dy;
        }

      }
    }

    if (turn_on && !moved.empty()) {
      handler->construct();
    }
    // Set the flag.
    turn_on = false;
  }

  void DemonWall::turnOn() {
    turn_on = true;
  }

}