#ifndef __ALIGNED_ARRAY_HPP__
#define __ALIGNED_ARRAY_HPP__

namespace GFlow {
  
  /*
   * @class aligned_array
   * Stores aligned data, like a vector
   */
  template<typename T> class aligned_array {
  public:
    // Default constructor
    aligned_array();
    
    // Presized array constructor
    aligned_array(int);
    
    // Fill array constructor
    aligned_array(int, const T&);
    
    // Destructor
    ~aligned_array();

    // Resize the array
    void reserve(int);

    // Resize the array in two parts
    void reserve2(int, int, int, int);
    
    // Data access
    T& at(int);
    T  at(int) const;
    T& operator[] (int);
    T  operator[] (int) const;
    
    // Get pointer
    T* getPtr() { return data; }

    // Accessors
    int size() { return _size; }

    // Mutators
    void setAlignment(int);

    // Exception classes: Out of bounds
    class aligned_array_out_of_bounds {
    public:
      aligned_array_out_of_bounds(int i) : index(i) {};
      int index;
    };

    // Exception classes: bad choice of alignment
    class bad_alignment {
    public:
      bad_alignment(int a) : align(a) {};
      int align;
    };
      
  private:
    int _size;
    int _alignment;
    
    T* data;
  };

#include "aligned_array.cpp"
  
}
#endif // __ALIGNED_ARRAY_HPP__
