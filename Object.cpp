#include "Object.h"

double Particle::getRatio() { 
  if (acceleration*acceleration == 0 || velocity*velocity == 0) return 1.0;
  return sqrt((velocity*velocity)/(acceleration*acceleration)); 
}

void Sphere::interact(Particle* P) {
  // Exert this particles force on both particles (in opposite directions of course)

  vect<> displacement = P->getPosition() - position;
  double distSqr = sqr(displacement);
  double cutoff = radius + P->getRadius();
  double cutoffsqr = sqr(cutoff);

  if (distSqr < cutoffsqr) { // Interaction
    double dist = sqrt(distSqr);
    vect<> normal = (1.0/dist) * displacement;

    // Repulsive force
    vect<> relV = P->getVelocity() - velocity;
    double strength = repulsion*(cutoff-dist)/P->getRadius();

    vect<> F = strength*normal;

    //** Got rid of dissipative forces for now
    //if (relV*displacement<0) F -= dissipation*(relV*normal)*normal;

    P->accelerate(F);
  }
}

void Sphere::update(double epsilon) {
  // Drag force
  vect<> dragF = -drag*velocity*velocity*normalize(velocity);
  accelerate(dragF);

  // Update velocity and position (velocity verlet)
  position += epsilon*(velocity + 0.5*epsilon*acceleration);
  velocity += 0.5*epsilon*(acceleration+acceleration_t);

  acceleration_t = acceleration; // Save last acceleration
  acceleration = vect<>(); // Reset acceleration
}

Wall::Wall(vect<> origin, vect<> wall) : origin(origin), wall(wall), coeff(0), repulsion(wall_repulsion), dissipation(wall_dissipation) {
  normal = wall;
  normal.normalize();
  length = wall.norm();
}

Wall::Wall(vect<> origin, vect<> end, bool) : origin(origin), wall(end-origin), coeff(0), repulsion(wall_repulsion), dissipation(wall_dissipation) {
  normal = wall;
  normal.normalize();
  length = wall.norm();
}

void Wall::interact(Particle* P) {

  vect<> displacement = P->getPosition() - origin;
  double l_par = displacement*normal;
  vect<> d_par = l_par*normal;
  vect<> d_perp = displacement - d_par;
  
  // Check whether the particle is between the start and end of the wall
  double radius = P->getRadius();
  double radSqr = sqr(radius);
  
  if (l_par>=0) { // Located forward of the origin
    if (length>l_par) displacement = d_perp;  // The particle is above the wall (in the perp. direction)
    else displacement -= wall; // Displacement from the nearest end (the far end) of the wall    
  }

  double distSqr = sqr(displacement);   // Located behind the origin

  /// We now have the correct displacement vector and distSqr value
  if (distSqr<=radSqr) {
    double dist = sqrt(distSqr);
    vect<> norm = normalize(displacement); // <normal> is already taken

    vect<> relV = P->getVelocity();
    double strength = repulsion*(radius-dist)/radius;

    vect<> F = strength*norm;
    if (relV*displacement<0) F -= dissipation*(norm*relV)*norm; // Only dissipate in the normal direction
    
    P->accelerate(F);
  }
}
