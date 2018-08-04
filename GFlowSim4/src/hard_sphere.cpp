#include "hard_sphere.hpp"
// Other files
#include "simdata.hpp"
#include "vectormath.hpp"

namespace GFlowSimulation {

  HardSphere::HardSphere(GFlow *gflow) : Force(gflow), repulsion(DEFAULT_HARD_SPHERE_REPULSION) {
    /*
    // --- TESTS
    BLOCK_SIZE = 128;
    const int n_particles = BLOCK_SIZE / DIMENSIONS; // Rounds down
    _x1   = new RealType[BLOCK_SIZE];
    _x2   = new RealType[BLOCK_SIZE];
    _disp = new RealType[BLOCK_SIZE];
    _dist = new RealType[n_particles];
    _f    = new RealType[BLOCK_SIZE];
    */
  };

  HardSphere::~HardSphere() {
    /*
    // --- TESTS
    if (_x1)   delete [] _x1;
    if (_x2)   delete [] _x2;
    if (_disp) delete [] _disp;
    if (_dist) delete [] _dist;
    if (_f)    delete [] _f;
    */
  }

  void HardSphere::calculateForces() const {
    int nverlet = verletList.vlSize(), id1(0), id2(0); // List length, id pointers
    if (nverlet==0) return; // No forces to calculate

    // Get the data we need
    RealType **x = Base::simData->x, **f = Base::simData->f;
    RealType *sg = Base::simData->sg;
    int *type = Base::simData->type;
    RealType displacement[DIMENSIONS], normal[DIMENSIONS]; // To calculate displacement, normal vector
    Bounds bounds = Base::gflow->getBounds(); // Simulation bounds
    BCFlag boundaryConditions[DIMENSIONS]; 
    copyVec(Base::gflow->getBCs(), boundaryConditions); // Keep a local copy of the wrap frags

    // Get verlet list data
    const int *verlet = verletList.getVerlet();
    RealType F[DIMENSIONS];

    // --- Go through all particles
    for (int i=0; i<nverlet; i+=2) {
      id1 = verlet[i];
      id2 = verlet[i+1];
      // Get the displacement between the particles
      getDisplacement(x[id1], x[id2], displacement, bounds, boundaryConditions);
      // Check if the particles should interact
      RealType dsqr = sqr(displacement);
      if (dsqr < sqr(sg[id1] + sg[id2])) {
        RealType distance = sqrt(dsqr);
        scalarMultVec(1./distance, displacement, normal);
        // Calculate force strength
        forceStrength(F, normal, distance, id1, id2);
      }
    }

    // --- START TESTS
    /*
    // How many particles will fit in each block - each particle 
    // needs [DIMENSIONS] space
    const int n_particles = BLOCK_SIZE / DIMENSIONS; // Rounds down

    const int block_steps = n_particles*DIMENSIONS;
    // How many complete blocks we need
    const int n_blocks = (nverlet * DIMENSIONS) / BLOCK_SIZE; // Rounds down
    // How much space is left over
    const int n_extra = nverlet - n_blocks*n_particles;
    
    cout << n_particles << " " << n_blocks << " " << n_extra << endl;

    // Copy particle positions to _x1, _x2
    for (int i=0; i<n_particles; ++i) {

      //***
      if (2*i+1>=nverlet) {
        cout << "Big 0.";
        exit(0);
      }
      //***

      copyVec(x[verlet[ 2*i ]], &_x1[i*DIMENSIONS]);
      copyVec(x[verlet[2*i+1]], &_x2[i*DIMENSIONS]);
    }
    // Find particle displacements - should be vectorized
    for (int b=0; b<block_steps; ++b) {
      //***
      if (b>=BLOCK_SIZE) {
        cout << "Big A.";
        exit(0);
      }
      //***
      _disp[b] = _x1[b] - _x2[b];
    }
    // Find particle distances
    for (int i=0; i<n_particles; ++i) {
      // Unroll this
      RealType distance = 0;
      for (int d=0; d<DIMENSIONS; ++d) {
        //***
        if (i*DIMENSIONS+d>=BLOCK_SIZE) {
          cout << "Big B.\n";
          cout << i << " " << n_particles << " " << i*DIMENSIONS << " " << i*DIMENSIONS+d << " " << BLOCK_SIZE << endl;
          exit(0);
        }
        //***
        distance += sqr(_disp[i*DIMENSIONS+d]);
      }
      // Calculate the distance factor. It is not just the distance.
      _dist[i] = distance<0.1 ? 0 : repulsion*(0.1 - sqrt(distance));
    } 
    // Calculate components of force - should be vectorized
    // _f    -> components of force 
    // _dist -> distance
    // _disp -> components of displacement
    for (int b=0; b<block_steps; ++b) {
      int id = b/DIMENSIONS;
      //***
      if (b>=BLOCK_SIZE) {
        cout << "Big C.";
        exit(0);
      }
      //***
      //***
      if (id>=n_particles) {
        cout << "Big D.";
        exit(0);
      }
      //***

      // We have set _dist[id] = _mask[b]*repulsion*(0.1 - _dist[id])
      _f[b] += _dist[id]*_disp[b];
    }

    // Add force components to where they should be
    */
  }

  void HardSphere::setRepulsion(RealType r) { 
    repulsion = r; 
  }

  inline void HardSphere::forceStrength(RealType *F, const RealType *normal, const RealType distance, const int id1, const int id2) const {
    // This should make sure that forces are zero if either object is of type -1. This does not seem to add much (any?) overhead
    RealType c1 = Base::simData->type[id1]<0 ? 0 : 1.; //--
    RealType c2 = Base::simData->type[id2]<0 ? 0 : 1.; //--

    scalarMultVec(c1*c2*repulsion*(simData->Sg(id1) + simData->Sg(id2) - distance), normal, F);
    // Add forces
    plusEqVec (Base::simData->f[id1], F);
    minusEqVec(Base::simData->f[id2], F);
  }

}
