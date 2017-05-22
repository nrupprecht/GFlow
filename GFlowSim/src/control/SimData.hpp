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

  // Forward declaration to DataRecord
  class DataRecord;
  // Forward declaration to ExternalForce
  class ExternalForce;
  // Forward declaration to Creator
  class Creator;
  // Forward declaration to 
  class FileParser;
  
  /*
   * @class SimData
   *
   */
  class SimData {
  public:
    // Constructor taking domain and simulation dimensions
    SimData(const Bounds&, const Bounds&);

    // Destructor
    ~SimData();

    // Reserve particle size
    void reserve(int, int=-1);

    // Reserve additional space for particles
    void reserveAdditional(int, int=-1);

    // Add walls
    void addWall(const Wall&);
    void addWall(const Bounds&);

    // Add particles
    void addParticle(const Particle&);
    void addParticle(const vector<Particle>&);

    // Access sizes and capacities
    int getDomainSize()     { return domain_size; }
    int getDomainCapacity() { return domain_capacity; }
    int getEdgeSize()       { return edge_size; }
    int getEdgeCapacity()   { return edge_capacity; }

    // Returns the bounds
    Bounds getBounds()    { return bounds; }
    Bounds getSimBounds() { return simBounds; } 

    // Return wrapping
    bool getWrapX()      { return wrapX; }
    bool getWrapY()      { return wrapY; }

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
    RealType* getRpPtr() { return rp.getPtr(); }
    RealType* getDsPtr() { return ds.getPtr(); }
    RealType* getCfPtr() { return cf.getPtr(); }
    int*      getItPtr() { return it.getPtr(); }

    // Summary pointer
    RealType** getPData();

    // Array access (non-checking)
    RealType& getPx(int i) { return px[i]; }
    RealType& getPy(int i) { return py[i]; }
    RealType& getVx(int i) { return vx[i]; }
    RealType& getVy(int i) { return vy[i]; }
    RealType& getSg(int i) { return sg[i]; }
    int&      getIt(int i) { return it[i]; }

    // Get position record
    vec2* getPRPtr()     { return positionRecord.getPtr(); }

    // Get walls
    vector<Wall>& getWalls() { return walls; }

    // Get a list of all the particles
    vector<Particle> getParticles();

    // Get the external forces
    vector<ExternalForce*>& getExternalForces() { return externalForces; }

    // Wrap positions and angles
    void wrap(RealType&, RealType&); // Position wrapping
    void wrap(RealType&);            // Angle wrapping

    // Displacement functions
    vec2 getDisplacement(const RealType, const RealType, const RealType, const RealType);
    vec2 getDisplacement(const vec2, const vec2);

#if USE_MPI == 1 // Get MPI data
    int getRank()    { return rank; }
    int getNumProc() { return numProc; }

    void atomMove();
    void atomCopy();
#endif

    /** Mutators **/

    // Set wrapping boundary conditions
    void setWrap(bool b)  { wrapX = wrapY = b; }
    void setWrapX(bool b) { wrapX = b; }
    void setWrapY(bool b) { wrapY = b; }

    // Set bounds
    void setBounds(const Bounds& b)    { bounds = b; }
    void setSimBounds(const Bounds& b) { simBounds = b; }

    // Add an external force or remove all of them
    void addExternalForce(ExternalForce*);
    void clearExternalForces();

    // DataRecord is a friend class
    friend class DataRecord;
    // Creator is a friend class
    friend class Creator;
    // FileParser is a friend class
    friend class FileParser;
    
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
    aligned_array<RealType> rp;
    aligned_array<RealType> ds;
    aligned_array<RealType> cf;
    // Interaction is -1 for empty array spaces
    aligned_array<int> it;

    // Data summary pointer
    RealType *pdata[15];
    void setPData();

    // Positions at the last verlet list creation
    aligned_array<vec2> positionRecord;

    // Wall data
    vector<Wall> walls;

    // Domain and Simulation bounds
    Bounds bounds;
    Bounds simBounds;

    // External forces
    vector<ExternalForce*> externalForces;

    // Simulation boundary conditions
    bool wrapX, wrapY;

#ifdef USE_MPI // Rank and number of processors - for MPI
    int rank, numProc;
#endif
    
  };

}
#endif // __SIM_DATA_HPP__
