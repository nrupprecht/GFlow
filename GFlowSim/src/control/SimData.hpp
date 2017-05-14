/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __SIM_DATA_HPP__
#define __SIM_DATA_HPP__

// Includes
#include "../../include/Utility.hpp"
#include "../../include/Bounds.hpp"
#include "../../include/aligned_array.hpp"
#include "../objects/Particle.hpp"
#include "../objects/Wall.hpp"

namespace GFlow {

  /*
   * @class SimData
   *
   */
  class SimData {
  public:
    // Constructor taking domain and simulation dimensions
    SimData(const Bounds&, const Bounds&);

    // Reserve particle size
    void reserve(int, int);

    // Add walls
    void addWall(const Wall&);

    // Add particles
    void addParticle(const Particle&);

    // Access sizes and capacities
    int getDomainSize()     { return domain_size; }
    int getDomainCapacity() { return domain_capacity; }
    int getEdgeSize()       { return edge_size; }
    int getEdgeCapacity()   { return edge_capacity; }

    // Returns the bounds
    Bounds getBounds()    { return bounds; }
    Bounds getSimBounds() { return simBounds; } 

    // Pointer access to arrays
    RealType* getPxPtr() { return px.getPtr(); }
    RealType* getPyPtr() { return py.getPtr(); }
    RealType* getVxPtr() { return vx.getPtr(); }
    RealType* getVyPtr() { return vy.getPtr(); }
    RealType* getFxPtr() { return fx.getPtr(); }
    RealType* getFyPtr() { return fy.getPtr(); }
    RealType* getThPtr() { return th.getPtr(); }
    RealType* getOmPtr() { return om.getPtr(); }
    RealType* getSgPtr() { return sg.getPtr(); }
    RealType* getTqPtr() { return tq.getPtr(); }
    RealType* getImPtr() { return im.getPtr(); }
    RealType* getIiPtr() { return iI.getPtr(); }
    int* getItPtr()      { return it.getPtr(); }

    // Get walls
    vector<Wall> getWalls() const { return walls; }

    // Wrap positions and angles
    void wrap(RealType&, RealType&); // Position wrapping
    void wrap(RealType&);            // Angle wrapping

    // Get MPI data
    int getRank()    { return rank; }
    int getNumProc() { return numProc; }

#ifdef USE_MPI
    void atomMove();
    void atomCopy();
#endif

  private:

    // Number of domain particles
    int domain_size;

    // Storage allocated for main particles
    int domain_capacity;

    // Number of edge particles
    int edge_size;
    
    // Storage allocated for edge particles
    int edge_capacity;

    // Particle data
    aligned_array<RealType> px;
    aligned_array<RealType> py;
    aligned_array<RealType> vx;
    aligned_array<RealType> vy;
    aligned_array<RealType> fx;
    aligned_array<RealType> fy;
    aligned_array<RealType> th;
    aligned_array<RealType> om;
    aligned_array<RealType> tq;
    aligned_array<RealType> sg;
    aligned_array<RealType> im;
    aligned_array<RealType> iI;
    // ... Other data ...
    aligned_array<int> it;

    // Wall data
    vector<Wall> walls;

    // Domain and Simulation bounds
    Bounds bounds;
    Bounds simBounds;

    // Simulation boundary conditions
    bool wrapX, wrapY;

    // Rank and number of processors - for MPI
    int rank, numProc;
    
  };

}
#endif // __SIM_DATA_HPP__
