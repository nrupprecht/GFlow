// Header file for Contiguous List (CL) data structure
// A list-like class. You can remove arbitrary entries in constant time and add entries to the end in vector-like time. All data entries are stored in a contiguous array to help with cache efficiency
// First Created: Feb 27, 2017

#ifndef CONTIGUOUSLIST_H
#define CONTIGUOUSLIST_H

// Unnamed namespace
namespace {
  const int default_list_length = 5;

  // List node structure for CL
  template <typename T> struct ContiguousListNode {
    ContiguousListNode() : data(T()), prev(-1), next(-1) {};
    T data;
    int prev, next; // Values of -1 means no predescesor/successor, -2 means erased
  };
}

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
  T& operator[] (int);

  class iterator {
  public:
    // Constructor
    iterator(ContiguousList* l) { array=l->array; point=0; lst=l; }
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
      if (next!=-1) array[next].prev = prev;
      array[point].data.~T(); // Call destructor
      array[point].next = array[point].prev = -2;
      lst->_size--;
      lst->_removals++;
    }

  private:
    int point;
    ContiguousListNode<T> *array;
    ContiguousList *lst; // Pointer to the list
  };

  class const_iterator {
  public:
    // Constructor
    const_iterator(ContiguousList* l) { array=l->array; point=0; lst=l; }
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
    // Equals operator
    const_iterator& operator= (const iterator& it) { point=it.point; array=it.array; }
    // Check whether this is the CL we are an iterator for
    bool check(ContiguousListNode<T> *a) const { return a==array; }
  private:
    int point;
    ContiguousListNode<T> *array;
    ContiguousList<T> *lst; // Pointer to the list
  };

  iterator begin() { return iterator(this); }
  const_iterator begin() const { return const_iterator(this); }
  iterator end() { return iterator(this, -1); }
  const_iterator end() const { return const_iterator(this, -1); }

  /// Mutators
  // Add an element to the list
  void push_back(const T& data);
  // Remove an entry by iterator
  void erase(iterator);

 private:
  ContiguousListNode<T> *array; // The stored data
  int _begin, _end;  // The first and last valid entries of array
  int _size;   // How many members the contiguous list has
  int _length; // How long array is
  
  int _removals; // How many elements have been removed (i.e. how many holes there are)
};

#include "ContiguousList.cpp"

#endif
