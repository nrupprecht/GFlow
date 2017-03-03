// Header file for Contiguous List (CL) data structure
// A list-like class. You can remove arbitrary entries in constant time and add entries to the end in vector-like time. All data entries are stored in a contiguous array to help with cache efficiency
// First Created: Feb 27, 2017

#ifndef CONTIGUOUSLIST_H
#define CONTIGUOUSLIST_H

// Unnamed namespace
//namespace {
  const int default_list_length = 5;

  // List node structure for CL
  template <typename T> struct ContiguousListNode {
    ContiguousListNode() : data(T()), prev(-1), next(-1) {};
    T data;
    int prev, next; // Values of -1 means no predescesor/successor, -2 means erased
  };
//}

// A list-like class. You can remove arbitrary entries in constant time and add entries to the end in vector-like time. All data entries are stored in a contiguous array to help with cache efficiency
template <typename T> class ContiguousList {
 public:
  // Constructor
  ContiguousList();
  // Set aside some number of entries in the array (initial length) 
  ContiguousList(int);
  
  /// Error classes
  class BadAccess {};    // Out of array bounds
  class DeletedEntry {}; // Trying to access a deleted entry

  /// Accessors
  int size() const { return _size; }
  int length() const { return _length; }
  int removals() const { return _removals; }
  bool empty() const { return _size==0; }
  T& operator[] (int);

  class iterator {
  public:
    // Constructor
    iterator(ContiguousList* l) { array=l->array; point=l->_begin; lst=l; }
    // Constructor
    iterator(ContiguousList* l, int pnt) { array=l->array; point=pnt; lst=l; }
    // Dereferencing operator
    T& operator* () { return array[point].data; }
    // Arrow operator
    T& operator-> () { return array[point].data; }
    // Preincrementation
    iterator operator++ () { int i=point; point=array[point].next; return iterator(lst, i); }
    // Postincrementation
    iterator operator++ (int) { point=array[point].next; return *this; }
    // Bool Equals operator
    bool operator== (const iterator& it) const { return array==it.array && point==it.point; }
    // Not equals
    bool operator!= (const iterator& it) const { return array!=it.array || point!=it.point; }
    // Equals operator 
    iterator& operator= (const iterator& it) { point=it.point; array=it.array; }
    // Check whether this is the CL we are an iterator for
    bool check(ContiguousListNode<T> *a) const { return a==array; }
    // Erase the element we are pointing at
    void erase() {
      int prev = array[point].prev, next = array[point].next;
      if (prev!=-1) array[prev].next = next;
      else _begin=next; // We erased the head node
      if (next!=-1) array[next].prev = prev;
      else _end=prev; // We erased the tail node
      array[point].data.~T(); // Call destructor
      array[point].next = array[point].prev = -2;
      lst->_size--;
      lst->_removals++;
    }
    // Pointer access
    T* getAddress() { return &array[point].data; }
    
  private:
    int point;
    ContiguousListNode<T> *array;
    ContiguousList *lst; // Pointer to the list
  };

  class const_iterator {
  public:
    // Constructor
    const_iterator(ContiguousList* l) { array=l->array; point=l->_begin; lst=l; }
    // Constructor
    const_iterator(ContiguousList* l, int pnt) { array=l->array; point=pnt; lst=l; }
    // Dereferencing operator
    T operator* () { return array[point].data; }
    // Arrow operator
    T operator-> () { return array[point].data; }
    // Preincrementation
    const_iterator operator++ () { int i=point; point=array[point].next; return iterator(lst, i); }
    // Postincrementation
    const_iterator operator++ (int) { point=array[point].next; return *this; }
    // Bool Equals operator
    bool operator== (const iterator& it) const { return array==it.array && point==it.point; }
    // Not equals
    bool operator!= (const iterator& it) const { return array!=it.array || point!=it.point; }
    // Check whether this is the CL we are an iterator for
    bool check(ContiguousListNode<T> *a) const { return a==array; }
    // Pointer access
    T* getAddress() { return &array[point].data;}

  private:
    int point;
    ContiguousListNode<T> *array;
    ContiguousList<T> *lst; // Pointer to the list
  };

  iterator begin() { return iterator(this); }
  iterator rbegin() { return iterator(this, _end); }
  const_iterator begin() const { return const_iterator(this); }
  const_iterator rbegin() const { return const_iterator(this, _end); }
  iterator end() { return iterator(this, -1); }
  const_iterator end() const { return const_iterator(this, -1); }

  /// Mutators
  // Add an element to the list. True if a resize occured
  bool push_back(const T& data);
  // Remove an entry by iterator
  void erase(iterator);
  // Clear the list
  void clear();
  // Delete the entry containing the T object pointed to by a T pointer
  void remove(T*);

 private:
  ContiguousListNode<T> *array; // The stored data
  int _begin, _end;  // The first and last valid entries of array
  int _size;   // How many members the contiguous list has
  int _length; // How long array is
  
  int _removals; // How many elements have been removed (i.e. how many holes there are)
};

#include "ContiguousList.cpp"

#endif
