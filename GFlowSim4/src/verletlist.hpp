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

    // Add a pair of interacting particles
    void addPair(const int, const int);

    // Set sizes (but not capacities) to zero, effectively "clearing" out the data
    void clear();

    // Return the total length of the verlet list 
    int vlSize() const;

    // Get a (const) pointer to the verlet array
    const int* getVerlet() const;

  private:

    // --- Helper functions
    inline void resizeVerlet(); 

    // --- Data
    int *verlet;
    int vsize, vcapacity;

    // --- For iteration through pairs of interacting particles
    mutable int _current_point;    // The current address in [verlet] that we are to look at
  };

}
#endif // __VERLET_LIST_HPP__GFLOW__
