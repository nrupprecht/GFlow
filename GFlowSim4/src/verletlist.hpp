#ifndef __VERLET_LIST_HPP__GFLOW__
#define __VERLET_LIST_HPP__GFLOW__ 

namespace GFlowSimulation {

  class VerletList {
  public:
    // Constructor
    VerletList();

    // Destructor
    ~VerletList();

    // Add a new head
    void addHead(int);

    // Add an element to the head
    void addToHead(int);

    // Set sizes (but not capacities) to zero, effectively "clearing" out the data
    void clear();

     // Return the last head added to the head array
    int lastHead();

  private:
    // --- Helper functions
    inline void resizeVerlet(); 
    inline void resizeHeads();

    // --- Data
    int *verlet, *heads;
    int vsize, hsize, vcapacity, hcapacity;
    int default_verlet_capacity, default_head_capacity;
  };

}
#endif // __VERLET_LIST_HPP__GFLOW__