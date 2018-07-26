#ifndef __VERLET_LIST_HPP__GFLOW__
#define __VERLET_LIST_HPP__GFLOW__ 

#include <functional>

namespace GFlowSimulation {

  class VerletList {
  public:
    // Constructor
    VerletList();

    // Copy Constructor
    VerletList(const VerletList&);

    // Destructor
    ~VerletList();

    // Equals operator
    VerletList& operator=(const VerletList&);

    // --- Mutators

    // Add a new head
    void addHead(int);

    // Add an element to the head
    void addToHead(int);

    // Set sizes (but not capacities) to zero, effectively "clearing" out the data
    void clear();

    // --- Accessors
    //! Sets up the verlet list for iteration through the pairs of particles. Sets the argument to be the first head.
    //! Returns false if the verlet list is empty.
    bool begin(int&) const;

    //! Allows us to iterate through pairs of interacting particles, sets the arguments to the head and second particle.
    //!  Returns false if we have reached the end of the verlet list
    bool next(int&, int&) const;

    //! Loops through pairs of interacting particles, passing them to the force interaction kernel
    void forceLoop(const class Force*) const;

    // Return the last head added to the head array
    int lastHead() const;

    // Return the total length of the verlet list 
    int vlSize() const;

    // Return the number of heads in the verlet list
    int vlHSize() const;

    // Get a (const) pointer to the verlet array
    const int* getVerlet() const;

    // Get a (const) pointer to the heads array
    const int* getHeads() const;

  private:

    // --- Helper functions
    inline void resizeVerlet(); 
    inline void resizeHeads();

    // --- Data
    int *verlet, *heads;
    int vsize, hsize, vcapacity, hcapacity;

    // --- For iteration through pairs of interacting particles
    mutable int _current_point;    // The current address in [verlet] that we are to look at
  };

}
#endif // __VERLET_LIST_HPP__GFLOW__
