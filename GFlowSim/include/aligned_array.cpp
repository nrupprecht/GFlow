// To be included in aligned_array.h
  
template<typename T> aligned_array<T>::aligned_array() : _size(0), _length(0), _alignment(64), data(0) {};

template<typename T> aligned_array<T>::aligned_array(int s) : _size(0), _length(1.5*s), _alignment(64) {
  data = (T*)aligned_alloc(_alignment, _length*sizeof(T));
  // Initialize the entries to default values
  for (int i=0; i<s; ++i) data[i] = T();
}

template<typename T> aligned_array<T>::aligned_array(int s, const T& d) : _size(0), _length(1.5*s), _alignment(64) {
  data = (T*)aligned_alloc(_alignment, _length*sizeof(T));
  int i;
  // Set the first s entries to the requested value
  for (i=0; i<s; ++i)   data[i] = d;
  // Initialize the rest of the entries to default values
  for (; i<length; ++i) data[i] = T();
}

template<typename T> aligned_array<T>::~aligned_array() {
  if (data) free(data);
  data = 0;
}

template<typename T> void aligned_array<T>::reserve(int s) {
  T* new_data = (T*)aligned_alloc(_alignment, s*sizeof(T));
  if (data) {
    int end = size<s ? size : s;
    for (int i=0; i<end; ++i) new_data[i] = data[i];
    free(data);
  }
  data = new_data;
}

template<typename T> T& aligned_array<T>::at(int i) {
  if (i<0 || _size<=i) throw aligned_array_out_of_bounds(i);
  return data[i];
}

template<typename T> T aligned_array<T>::at(int i) const {
  if (i<0 || _size<=i) throw aligned_array_out_of_bounds(i);
  return data[i];
}

template<typename T> T& aligned_array<T>::operator() (int i) {
  return data[i];
}

template<typename T> T aligned_array<T>::operator() (int i) const {
  return data[i];
}

template<typename T> void aligned_array<T>::setAlignment(int a) {
  if (a<0) throw bad_alignment(a);
  _alignment = a;
  reserve(_length);
}
