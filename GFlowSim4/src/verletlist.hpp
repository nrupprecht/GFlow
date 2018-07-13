#ifndef __VERLET_LIST_HPP__GFLOW__
#define __VERLET_LIST_HPP__GFLOW__ 

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
    bool begin(int&);

    //! Allows us to iterate through pairs of interacting particles, sets the arguments to the head and second particle.
    //!  Returns false if we have reached the end of the verlet list
    bool next(int&, int&);

    //! Loops through pairs of interacting particles, passing them to the force interaction kernel
    void forceLoop(class Force*);

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
    int default_verlet_capacity, default_head_capacity;

    // --- For iteration through pairs of interacting particles
    int _current_point;    // The current address in [verlet] that the second particle (id2) is
    int _next_head;        // The address in [verlet] where the next head is
    int _next_head_number; // The address in [heads] where the next head is
    bool _last_region;     // True if we are iterating through the list for the last head particle
  };

}
#endif // __VERLET_LIST_HPP__GFLOW__