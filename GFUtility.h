#ifndef GFUTILITY_H
#define GFUTILITY_H

#include "Utility.h"
#include "Shape.h"

/// Array copy function
template<typename T> inline void copy(T *source, T *destination, int num) {
  for (int i=0; i<num; i++) destination[i] = source[i];
}

/// Array create and copy function
template<typename T> inline void create(T *source, T*& destination, int num) {
  destination = new T[num];
  copy(source, destination, num);
}

/// Arbitrary length storage
template<typename Type> class glStorage {
 public:
  glStorage() {
    total = 0;
    entries = 0;
  }
  glStorage(Type first) {
    total = 1;
    entries = new Type(first);
  }
  glStorage(const glStorage<Type>& store) {
    total = store.total;
    entries = new Type[total];
    for (int i=0; i<total; i++) entries[i] = store.entries[i];
  }
  glStorage(glStorage<Type>&& store) {
    entries = store.entries;
    total = store.total;
  }
  template<typename ...T> glStorage(int first, T... last) {
    vector<int> vec;
    total = sizeof...(last) + 1;
    entries = new Type[total];
    if (total==0) {
      entries = 0;
      return;
    }
    getEntries(0, first, last...);
  }

  glStorage& operator=(const glStorage<Type>& store) {
    if (store.total!=total) {
      total = store.total;
      if (entries)
	delete [] entries;
      entries = new Type[total];
    }
    for (int i=0; i<total; i++) entries[i] = store.entries[i];
  }

  glStorage& operator=(glStorage<Type>&& store) {
    total = store.total;
    if (entries) delete [] entries;
    entries = store.entries;
  }

  bool operator==(const glStorage<Type>& store) const {
    if (store.total!=total) return false;
    for (int i=0; i<total; i++)
      if (store.entries[i]!=entries[i]) return false;
    return true;
  }

  friend ostream& operator<<(ostream& out, glStorage<Type> store) {
    out << '{';
    for (int i=0; i<store.total; i++) {
      out << store.entries[i];
      if (i!=store.total-1) out << ',';
    }
    out << '}';
    return out;
  }

  glStorage operator+(const glStorage<Type>& store) const {
    if (store.total!=total) throw glStorageMismatch();
    glStorage<Type> in = store;
    for (int i=0; i<store.total; i++) in.entries[i] += entries[i];
    return in;
  }

  glStorage operator-(const glStorage<Type>& store) const {
    if (store.total!=total) throw glStorageMismatch();
    glStorage<Type> in = *this;
    for (int i=0; i<store.total; i++) in.entries[i] -= store.entries[i];
    return in;
  }

  int size() const { return total; }

  Type& at(int i) { return entries[i]; }
  Type at(int i) const { return entries[i]; }

  void resize(int s) {
    if (entries) {
      if (s!=total) {
	total = s;
	delete [] entries;
	entries = new Type[total];
	for (int i=0; i<total; i++) entries[i] = 0;
      }
    }
    else {
      total = s;
      entries = new Type[total];
      for (int i=0; i<total; i++) entries[i] = 0;
    }
  }

  class glStorageMismatch {};
  class glMismatch {};

 private:
  /// Helper functions
  template<typename ...T> void getEntries(int step, int first, T... last) const {
    entries[step] = first;
    getEntries(step+1, last...);
  }
  void getEntries(int step, int first) const {
    entries[step] = first;
  }

  /// Data
  Type *entries;
  int total;
};

template<typename T> inline string bare_representation(const glStorage<T>& gls) {
  stringstream stream;
  string str;
  for (int i=0; i<gls.size(); i++) {
    stream << gls.at(i);
    if (i!=gls.size()-1) stream << ',';
  }
  stream >> str;
  return str;
}

typedef glStorage<int> Index;
typedef glStorage<double> gVector;

typedef double (*gFunction) (gVector);

#endif
