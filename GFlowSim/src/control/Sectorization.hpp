/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 11, 2017
 *
 */

#ifndef __SECTORIZATION_HPP__
#define __SECTORIZATION_HPP__

// Includes
#include "../../include/Utility.hpp"
#include "SimData.hpp"

namespace GFlow {

  // Forward declaration to DataRecord
  class DataRecord;

  /*
   * @class Sectorization
   * Uses cell lists to create verlet lists so we can easily calculate forces
   */
  class Sectorization {
  public:
    // Default constructor
    Sectorization();

    // Sim data constructor
    Sectorization(SimData*);

    // Destructor
    ~Sectorization();
    
    // Set the simulation data to manage, set up sectorization
    void initialize(SimData* sd);
    
    // Update sectorization - calls the virtual function _sectorize()
    void sectorize();

    // Check if we need to remake our lists
    RealType checkNeedRemake();

    // Create verlet lists (bool is for force make list)
    void createVerletLists();

    // Create wall lists (bool is for force make list)
    void createWallLists();

    // Remove particles that overlap by to much
    void removeOverlapping(RealType=0.25);

    /****** Accessors *****/

    // Returns the particle verlet list
    auto& getVerletList() { return verletList; }

    // Returns the wall neighbor list
    auto& getWallList() { return wallList; }

    // Displacement functions
    vec2 getDisplacement(const RealType, const RealType, const RealType, const RealType);
    vec2 getDisplacement(const vec2, const vec2);

    // Data accessors
    int getNSX()            { return nsx; }
    int getNSY()            { return nsy; }
    RealType getSDX()       { return sdx; }
    RealType getSDY()       { return sdy; }
    RealType getCutoff()    { return cutoff; }
    RealType getSkinDepth() { return skinDepth; }
    RealType getMovement()  { return movement; }
    RealType getMaxCutR()   { return maxCutR; }
    RealType getSecCutR()   { return secCutR; }

    // Get the closest
    int getClosest(int, SimData*);

    // Get closest two
    pair<int, int> getClosestTwo(int, SimData*);

    // DataRecord is a friend class
    friend class DataRecord;

  protected:
    virtual void _sectorize();
    virtual void _createVerletLists();
    virtual void _createWallLists();
    virtual void _makeSectors();

    // Private helper functions
    inline int getSec(const RealType, const RealType);
    inline vector<int>& sec_at(int x, int y) { return sectors[nsx*y+x]; }

    // Pointer to the simulation data we manage
    SimData *simData;
    
    // The sectors: a list of the id's of particles that are in these sectors
    vector<int> *sectors;

    // Particle neighbor list. First int is the id of the main particle, the other ints are particles that are withing the cutoff region of the main particle.
    VListType verletList;

    // Wall neighbor list. First int is the id of the wall, the other ints are the id's of the particles within the cutoff region.
    WListType wallList;

    // Sectorization Data
    RealType cutoff;     // The cutoff radius
    RealType skinDepth;  // The skin depth
    RealType movement;   // The max possible movement of particles relative to one another since the last verlet list construction
    RealType maxCutR, secCutR; // First and second largest particle cutoffs
    int nsx, nsy;        // The number of sectors
    RealType sdx, sdy;   // Sector width and height
    RealType isdx, isdy; // The inverse of the sector widths and heights

    // Data from SimData
    Bounds bounds;       // Domain bounds
    Bounds simBounds;    // Simulation bounds
    bool wrapX, wrapY;   // Boundary conditions

    // Rank and number of processors - for MPI 
    int rank, numProc;
  };

}
#endif // __SECTORIZATION_HPP__
