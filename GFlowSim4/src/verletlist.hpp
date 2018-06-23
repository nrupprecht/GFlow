#ifndef __VERLET_LIST_HPP__GFLOW__
#define __VERLET_LIST_HPP__GFLOW__ 

namespace GFlowSimulation {

  class VerletList {
  public:
    // Constructor
    VerletList();

    // Destructor
    ~VerletList();

    // --- Mutators

    // Add a new head
    void addHead(int);

    // Add an element to the head
    void addToHead(int);

    // Set sizes (but not capacities) to zero, effectively "clearing" out the data
    void clear();

    // --- Accessors

     // Return the last head added to the head array
    int lastHead();

    // Return the total length of the verlet list 
    int vlSize();

    // Return the number of heads in the verlet list
    int vlHSize();

    // Get a (const) pointer to the verlet array
    const int* getVerlet();

    // Get a (const) pointer to the heads array
    const int* getHeads();

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