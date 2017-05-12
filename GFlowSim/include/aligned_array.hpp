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
    
    // Data acces
    T& at(int);
    T  at(int) const;
    T& operator() (int);
    T  operator() (int) const;
    
    // Accessors
    int size() { return _size; }
    int length() { return _length; }

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
    int _length;
    int _alignment;
    
    T* data;
  };

#include "aligned_array.cpp"
  
}
#endif // __ALIGNED_ARRAY_HPP__
