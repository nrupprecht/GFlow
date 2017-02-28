// Included in ContiguousList.h

template<typename T> ContiguousList<T>::ContiguousList() : array(new ContiguousListNode<T>[default_list_length]), _begin(0), _end(0), _size(0), _length(default_list_length), _removals(0) {};

template<typename T> ContiguousList<T>::ContiguousList(int reserve) : array(new ContiguousListNode<T>[reserve]), _begin(0), _end(0), _size(0), _length(reserve), _removals(0) {};

template<typename T> T& ContiguousList<T>::operator[] (int i) {
  if (length<=i || i<0) throw BadAccess();
  if (array[i].next==-2) throw DeletedEntry();
  return array[i].data;
}

template<typename T> void ContiguousList<T>::push_back(const T& data) {
  if (_length==0) { // A pathalogical case
    array = new ContiguousListNode<T>[default_list_length];
    _length = default_list_length;
  }
  if (_size==0) { // Initial entry or a pathalogical case
    _begin = 0;
    _end = 0;
    _size = 1;
    array[0].data = data;
    array[0].next = array[0].prev = -1;
    return;
  }
  else if (_size==_length) { // Create new array (if length==0, we arrive here)
    _length = 1.5*_size + 1; // length > 0
    ContiguousListNode<T>* newArray = new ContiguousListNode<T>[_length];
    int point = _begin, count = 0;
    while (point!=_end) {
      newArray[count].data = array[point].data;
      newArray[count].next = count+1;
      newArray[count].prev = count-1; // For first element, prev = -1
      point = array[point].next;
      ++count;
    }
    // point == end
    newArray[count].data = array[_end].data;
    newArray[count].next = -1;
    newArray[count].next = count-1;
    // Set begin and end
    _begin = 0; _end = count;
    delete [] array;
    array = newArray;
    _removals = 0;
  }
  // Add data to the end of the array
  array[_end].next = _end+1;
  _end++;
  array[_end].data = data;
  array[_end].next = -1; // For last element, next = -1
  array[_end].prev = _end-1;
  _size++;
}

template<typename T> void ContiguousList<T>::erase(iterator it) {
  if (it.check(array)) it.erase();
}
