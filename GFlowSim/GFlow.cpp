#include "GFlow.h"

using MPI::COMM_WORLD;

void GFlow::createSquare(int N, double radius) {
  Bounds bounds(0, 1, 0, 1);
  setBounds(bounds);
  // Set cutoff and skin depth
  cutoff = 2*radius;
  skinDepth = 0.25*radius;
  // Create particles on the root processor
  list<Particle> allParticles;
  if (rank==0) {
    for (int i=0; i<N; ++i) {
      vect<> pos(drand48(), drand48());
      double angle = 2*PI*drand48();
      vect<> v(cos(angle), sin(angle));
      Particle p(pos, radius);
      p.velocity = v;
      allParticles.push_back(p);
    }
  }
  COMM_WORLD.Barrier();
  // Take your own particles
  if (rank==0) {
    Bounds b = getBoundsForProc(0);
    vector<list<Particle>::iterator> remove;
    for (auto p=allParticles.begin(); p!=allParticles.end(); ++p)
      if (b.contains(p->position)) {
	remove.push_back(p);
	particles.push_back(*p);
      }
    for (auto &p : remove) allParticles.erase(p);
  }
  // Distribute particles to processes
  for (int proc=1; proc<ndx*ndy; ++proc) {
    if (rank==0) { // Send
      vector<Particle> domainParticles;
      Bounds b = getBoundsForProc(proc);
      vector<list<Particle>::iterator> remove;
      for (auto p=allParticles.begin(); p!=allParticles.end(); ++p)
	if (b.contains(p->position)) {
	  remove.push_back(p);
	  domainParticles.push_back(*p);
	}
      for (auto &p : remove) allParticles.erase(p);
      int size = domainParticles.size();
      Particle *buffer = new Particle[size];
      for (int j=0; j<size; ++j) buffer[j] = domainParticles[j];
      // Send the amount of data we are going to send to processor i
      COMM_WORLD.Send( &size, 1, MPI_INT, proc, 0);
      // Send the actual data to processor i
      COMM_WORLD.Send( buffer, size, PARTICLE, proc, 0);
    }
    else if (rank==proc) { // Recieve
      int size = 0, root = 0;
      // Recieve the amount of data we should expect
      COMM_WORLD.Recv( &size, 1, MPI_INT, root, 0);
      // Recieve the actual data
      Particle *buffer = new Particle[size];
      COMM_WORLD.Recv( buffer, size, PARTICLE, root, 0);
    }
    else; // Nothing
  }
  COMM_WORLD.Barrier();
  
}
