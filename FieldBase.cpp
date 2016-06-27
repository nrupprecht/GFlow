// .cpp file for template functions

template<typename T>
FieldBase<T>::FieldBase() {
  initialize();
}

template<typename T>
FieldBase<T>::FieldBase(int x, int y) {
  initialize();
  setDims(x,y);
}

template<typename T>
void FieldBase<T>::initialize() {
  left = bottom = 0;
  right = top = 1;
  dX = dY = 0;
  invDist = Zero; // invDist is meaningless
  array = 0;
  wrapX = true; wrapY = false;
  usesLocks = false;
  locks = 0;
  solveIterations = 500;
  tollerance = 0.01;
}

template<typename T>
void FieldBase<T>::setDims(int x, int y) {
  wrapX = true; wrapY = false;
  dX = x; dY = y;
  if (wrapX) invDist.x = (right-left)/dX;
  else invDist.x = (right-left)/(dX-1);
  if (wrapY) invDist.y = (top-bottom)/dY;
  else invDist.y = (top-bottom)/(dY-1);

  LFactor = 2*(sqr(invDist.x)+sqr(invDist.y));
  invLFactor = 1./LFactor;

  array = new T[dX*dY];
  for (int i=0; i<dX*dY; i++) array[i] = T(0);
  if (usesLocks) createLocks();
}

template<typename T>
void FieldBase<T>::setWrapX(bool w) {
  wrapX = w;
  if (w) invDist.x = (right-left)/dX;
  else invDist.x = (right-left)/(dX-1);
  LFactor = 2*(sqr(invDist.x)+sqr(invDist.y));
  invLFactor = 1./LFactor;
}

template<typename T>
void FieldBase<T>::setWrapY(bool w) {
  wrapY = w;
  if (w) invDist.y = (top-bottom)/dY;
  else (top-bottom)/(dY-1);
  LFactor = 2*(sqr(invDist.x)+sqr(invDist.y));
  invLFactor = 1./LFactor;
}

template<typename T>
void FieldBase<T>::setEdges(double x) {
  for (int i=0; i<dX; i++) {
    at(i,0) = x;
    at(i,dY-1) = x;
  }
  for (int j=1; j<dY-1; j++) {
    at(0,j) = x;
    at(dX-1,j) = x;
  }
}

template<typename T>
T& FieldBase<T>::at(int x, int y) {
  if (wrapX) {
    while (x<0) x+=dX;
    while (x>=dX) x-=dX;
  }
  if (wrapY) {
    while (y<0) y+=dY;
    while (y>=dY) y-=dY;
  }
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds();
  return array[y*dX+x];
}

template<typename T>
T FieldBase<T>::at(int x, int y) const {
  if (wrapX) {
    while (x<0) x+=dX;
    while (x>=dX) x-=dX;
  }
  if (wrapY) {
    while (y<0) y+=dY;
    while (y>=dY) y-=dY;
  }
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds();
  return array[y*dX+x];
}

template<typename T>
T& FieldBase<T>::operator()(int x, int y) {
  return at(x,y);
}

template<typename T>
string FieldBase<T>::print() const {
  stringstream stream;
  stream << '{';
  for (int y=dY-1; y>=0; y--) {
    stream << '{';
    for (int x=0; x<dX; x++) {
      stream << at(x,y);
      if (x!=dX-1) stream << ',';
    }
    stream << '}';
    if (y!=0) stream << ',';
  }
  stream << '}';

  string str;
  stream >> str;
  return str;
}

template<typename T>
void FieldBase<T>::lock(int x, int y, bool l) {
  if (usesLocks) {
    locks[x+dX*y] = l;
  }
}

template<typename T>
bool& FieldBase<T>::lockAt(int x, int y) {
  static bool L = false;
  if (wrapX) {
    while (x<0) x+=dX;
    while (x>=dX) x-=dX;
  }
  if (wrapY) {
    while (y<0) y+=dY;
    while (y>=dY) y-=dY;
  }
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds();
  if (!usesLocks) return L;
  return locks[x+dX*y];
}

template<typename T>
void FieldBase<T>::lockEdges(bool l) {
  for (int i=0; i<dX; i++) {
    lockAt(i,0) = l;
    lockAt(i,dY-1) = l;
  }
  for (int j=1; j<dY-1; j++) {
    lockAt(0,j) = l;
    lockAt(dX-1,j) = l;
  }
}

template<typename T>
FieldBase<T> FieldBase<T>::operator+=(T x) {
  for (int i=0; i<dX*dY; i++) array[i] += x;
  return *this;
}

template<typename T>
void FieldBase<T>::plusEq(FieldBase& field, double mult) {
  for (int i=0; i<dX*dY; i++) array[i] += mult*field.array[i];
}

template<typename T>
void FieldBase<T>::minusEq(FieldBase& field, double mult) {
  for (int i=0; i<dX*dY; i++) array[i] -= mult*field.array[i];
}

template<typename T>
T FieldBase<T>::DX(int x, int y) const {
  if (wrapX) {
    if (x==0) return 0.5*invDist.x*(at(1,y)-at(dX-1,y));
    else if (x==dX-1) return 0.5*invDist.x*(at(0,y)-at(dX-2,y));
    else return 0.5*invDist.x*(at(x+1,y)-at(x-1,y));
  }
  else {
    if (y==0) return invDist.x*(at(1,y)-at(0,y));
    else if (y==dY-1) return invDist.x*(at(dX-1,y)-at(dX-2,y));
    else return 0.5*invDist.x*(at(x+1,y)-at(x-1,y));
  }
}

template<typename T>
T FieldBase<T>::DY(int x, int y) const {
  if (wrapY) {
    if (y==0) return 0.5*invDist.y*(at(x,1)-at(x,dY-1));
    else if (y==dY-1) return 0.5*invDist.y*(at(x,0)-at(x,dY-2));
    else return 0.5*invDist.y*(at(x,y+1)-at(x,y-1));
  }
  else {
    if (y==0) return invDist.y*(at(x,1)-at(x,0));
    else if (y==dY-1) return invDist.y*(at(x,dY-1)-at(x,dY-2));
    else return 0.5*invDist.y*(at(x,y+1)-at(x,y-1));
  }
}

template<typename T>
vect<> FieldBase<T>::getPos(int x, int y) const {
  return vect<>(left + x*invDist.x, bottom + y*invDist.y);
}

template<typename T>
void FieldBase<T>::resetLocks() {
  if (usesLocks && locks) {
    for (int i=0; i<dX*dY; i++) locks[i] = false;
  }
}

template<typename T>
void FieldBase<T>::useLocks() {
  usesLocks = true;
  createLocks();
}

template<typename T>
void FieldBase<T>::createLocks() {

  cout << locks << endl;

  if (locks) delete [] locks;
  locks = new bool[dX*dY];
  for (int i=0; i<dX*dY; i++) locks[i] = false;
}

/// Successive Over-Relaxation (SOR) using checkerboard updating with no source
template<typename T>
void FieldBase<T>::SOR_solver() {
  // Calculate relaxation parameter
  double omega = 2.0/(1+PI/dX);
  double maxDelta = 1e9;
  for (int iter=0; iter<solveIterations && maxDelta>tollerance; iter++) {
    // Even squares
    maxDelta = 0;
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
        if ((x+y)%2==0 && !lockAt(x,y)) {
          auto value = (1-omega)*at(x,y) + omega*invLFactor*(sqr(invDist.x)*(at(x+1,y)+at(x-1,y)) + sqr(invDist.y)*(at(x,y+1)+at(x,y-1)));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
    // Odd squares
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
	if ((x+y)%2 && !lockAt(x,y)) {
          auto value = (1-omega)*at(x,y) + omega*invLFactor*(sqr(invDist.x)*(at(x+1,y)+at(x-1,y)) + sqr(invDist.y)*(at(x,y+1)+at(x,y-1)));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
  }
}

/// Successive Over-Relaxation (SOR) using checkerboard updating with source
template<typename T>
void FieldBase<T>::SOR_solver(FieldBase& source, double mult) {
  // Check that the source field has the same dimensions
  if (source.dX!=dX || source.dY!=dY) throw FieldMismatch();
  // Calculate relaxation parameter
  double omega = 2.0/(1+PI/dX);
  double maxDelta = 1e9;
  for (int iter=0; iter<solveIterations && maxDelta>tollerance; iter++) {
    // Even squares
    maxDelta = 0;
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
        if ((x+y)%2==0 && !lockAt(x,y)) {
          auto value = (1-omega)*at(x,y) + omega*invLFactor*(sqr(invDist.x)*(at(x+1,y)+at(x-1,y)) + sqr(invDist.y)*(at(x,y+1)+at(x,y-1)) - mult*source(x,y));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
    // Odd squares
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
        if ((x+y)%2 && !lockAt(x,y)) {
          auto value = (1-omega)*at(x,y) + omega*invLFactor*(sqr(invDist.x)*(at(x+1,y)+at(x-1,y)) + sqr(invDist.y)*(at(x,y+1)+at(x,y-1)) - mult*source(x,y));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
  }
}
