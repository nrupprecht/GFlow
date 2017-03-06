#include "Sectorization.h"

Sectorization::Sectorization() : particles(0) {
  // MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  // MPI_Comm_size( MPI_COMM_WORLD, &numProc );

  // MPI::Init();
  rank = MPI::COMM_WORLD.Get_rank();
  numProc = MPI::COMM_WORLD.Get_size();  

  MPI_Datatype PARTICLE;
  MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
  MPI_Type_commit( &PARTICLE );

  //**
  // cout << rank << ", " << numProc << "\n";
}

void Sectorization::interactions() {
  for (int y=1; y<nsy-1; ++y) 
    for (int x=1 ; x<nsx-1 ; ++x) 
      for (auto &p : sectors[y*nsx+x]) {
	// Check for interactions with particles in the surrounding sectors
	
      }
}

void Sectorization::update() {
  double dt = 0.5 * epsilon;
  for (auto &p : particles) {
    double mass = 1./p.invMass;
    p.velocity += dt * p.invMass * p.force;
    p.position += epsilon * p.velocity;
    // Apply gravity (part of step 3)
    p.force += gravity*mass;
  }
  // Update sectorization, keep track of particles that need to migrate to other processors, send particles to the correct processor
  updateSectors();
  // Reset particle's force recordings
  for (auto &p : particles) p.force = Zero;
  // Interaction forces
  interactions();
  // Do behaviors (particle characteristics)
  
  // Velocity update part two (step four)
  for (auto &p :particles) p.velocity += dt * p.invMass * p.force;
}

void Sectorization::updateSectors() {
  
  int data = rank, recv = -1;
  // Even Send
  if (rank%2==0) {
    if (rank+1<numProc) {
      int destination = rank+1;
      MPI_Send( &data, 1,MPI_INT, destination, 0, MPI_COMM_WORLD);
    }
  }
  else {
    int sender = rank-1;
    MPI_Recv(&recv, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
  // Odd send
  if (rank%2) {
    int destination = (rank+1)%numProc;
    MPI_Send( &data, 1, MPI_INT, destination, 0, MPI_COMM_WORLD );
  }
  else {
    if (rank!=0 || numProc%2==0) {
      int sender = rank!=0 ? rank-1 : numProc-1;
      MPI_Recv( &recv, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }
  }
  // Wrap if numProc is odd
  if (numProc%2) {
    if (rank==numProc-1) {
      int destination = 0;
      MPI_Send( &data, 1, MPI_INT, destination, 0, MPI_COMM_WORLD);
    }
    else if (rank==0) {
      int sender = numProc-1;
      MPI_Recv( &recv, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }
  }
  
  cout << "Rank, recieved: " << rank << " " << recv << endl;
}

void Sectorization::updateVelocities() {
  
}

void Sectorization::updatePositions() {

}
