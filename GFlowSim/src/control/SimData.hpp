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
#include "../objects/Characteristic.hpp"

namespace GFlow {

  // Forward declaration to DataRecord
  class DataRecord;
  // Forward declaration to ExternalForce
  class ExternalForce;
  // Forward declaration to Creator
  class Creator;
  // Forward declaration to FileParser
  class FileParser;
  // Forward declaration to Sectorization
  class Sectorization;
  // Forward declaration to ForceHandler
  class ForceHandler;

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
    void addWall(const vector<Wall>&);

    // Add particles
    int addParticle(const Particle&);
    void addParticle(const vector<Particle>&);
    int addParticle(const Particle&, Characteristic*);

    // Add a characteristic to a particle
    void addCharacteristic(int, Characteristic*);

    // Remove particles
    void removeAt(int);

    // Access sizes and capacities
    int getDomainSize()     { return domain_size; }
    int getDomainEnd()      { return domain_end; }
    int getDomainCapacity() { return domain_capacity; }
    int getEdgeSize()       { return edge_size; }
    int getEdgeCapacity()   { return edge_capacity; }
    
    // Get the current time
    RealType getTime() { return time; }

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

    // Array access (non-checking)
    RealType& getPx(int i) { return px[i]; }
    RealType& getPy(int i) { return py[i]; }
    RealType& getVx(int i) { return vx[i]; }
    RealType& getVy(int i) { return vy[i]; }
    RealType& getFx(int i) { return fx[i]; }
    RealType& getFy(int i) { return fy[i]; }
    RealType& getOm(int i) { return om[i]; }
    RealType& getSg(int i) { return sg[i]; }
    RealType& getIm(int i) { return im[i]; }
    int&      getIt(int i) { return it[i]; }

    // Get a copy of a particle
    Particle makeParticle(int);

    // Get position record
    vec2* getPRPtr()     { return positionRecord.getPtr(); }

    // Get pressure data
    void getPressureData(vector<PData>&, RealType=0.);

    // Get walls
    vector<Wall>& getWalls() { return walls; }

    // Get a list of all the particles
    vector<Particle> getParticles();

    // Get all the particles withing a cutoff radius
    vector<Particle> getParticles(vec2, RealType);
    vector<int> getParticlesID(vec2, RealType);

    // Get the characteristics
    auto& getCharacteristics() { return characteristics; }

    // Get the force handler
    ForceHandler* getForceHandler() { return forceHandler; }

    // Get the external forces
    vector<ExternalForce*>& getExternalForces() { return externalForces; }

    // Get termination indicator
    bool getTerminate() { return terminate; }

    // Wrap positions and angles
    void wrap(RealType&, RealType&); // Position wrapping
    void wrap(RealType&);            // Angle wrapping

    // Displacement functions
    vec2 getDisplacement(const RealType, const RealType, const RealType, const RealType);
    vec2 getDisplacement(const vec2, const vec2);
    vec2 getWallDisplacement(const Wall&, const vec2, RealType);

    // Calculate packing fraction
    RealType getPhi();

    // Get the indices of the two closest particles to a particular particle
    pair<int, int> getClosestTwo(int);

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

    // Set sectorization
    void setSectors(Sectorization *sec) { sectors = sec; }

    // Set force handler
    void setForceHandler(ForceHandler* frc) { forceHandler = frc; }

    // Add an external force or remove all of them
    void addExternalForce(ExternalForce*);
    void clearExternalForces();

    // Clear forces and torques
    void clearForceTorque();

    // Update position record
    void updatePositionRecord();

    // Set initial position vector
    void setInitialPositions();

    // Get initial position vector
    const vector<vec2>& getInitialPositions() { return initialPositions; }

    // Remove particles that overlap too much
    void removeOverlapping(RealType=0.03);

    // Set termination indicator
    void setTerminate(bool s) { terminate = s; }

    // DataRecord is a friend class
    friend class DataRecord;
    // Creator is a friend class
    friend class Creator;
    // FileParser is a friend class
    friend class FileParser;
    // Integrator is a friend class
    friend class Integrator;

    struct BadParticle {
      BadParticle(int i) : index(i) {};
      int index;
    };
    
  private:
    // Private helper functions
    inline void compressArrays();

    // Current time
    RealType time;

    // Number of domain particles
    int domain_size;

    // Position after the last valid particle in the domain
    int domain_end;

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

    // Positions at the last verlet list creation
    aligned_array<vec2> positionRecord;

    // Wall data
    vector<Wall> walls;

    // Holes in the arrays
    list<int> holes;

    // Initial positions
    vector<vec2> initialPositions;

    // Characteristics - we are in charge of managing (deleting) these
    std::map<int, Characteristic*> characteristics;

    // Domain and Simulation bounds
    Bounds bounds;
    Bounds simBounds;

    // External forces
    vector<ExternalForce*> externalForces;

    // Simulation boundary conditions
    bool wrapX, wrapY;

    // The sectorization that handles this data
    Sectorization *sectors;

    // The force handler that handles our forces
    ForceHandler *forceHandler;

    // A termination indicator
    bool terminate;

#ifdef USE_MPI // Rank and number of processors - for MPI
    int rank, numProc;
#endif
    
  };

}
#endif // __SIM_DATA_HPP__
