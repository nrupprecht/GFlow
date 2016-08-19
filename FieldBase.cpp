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
  wrapX = true; wrapY = true;
  usesLocks = false;
  locks = 0;
  solveIterations = 500;
  tollerance = 0.0001;
}

template<typename T>
FieldBase<T>& FieldBase<T>::operator=(T val) {
  for (int y=0; y<dY; y++)
    for (int x=0; x<dX; x++)
      at(x,y) = val;
  return *this;
}

template<typename T>
void FieldBase<T>::setDims(int x, int y) {
  dX = x; dY = y;
  // Recompute distances
  if (wrapX) invDist.x = (right-left)/dX;
  else invDist.x = (right-left)/(dX-1);
  if (wrapY) invDist.y = (top-bottom)/dY;
  else invDist.y = (top-bottom)/(dY-1);
  LFactor = 2*(sqr(invDist.x)+sqr(invDist.y));
  invLFactor = 1./LFactor;
  // Create new array
  array = new T[dX*dY];
  for (int i=0; i<dX*dY; i++) array[i] = T(0);
  if (usesLocks) createLocks();
}

template<typename T>
void FieldBase<T>::setBounds(double l, double r, double b, double t) {
  left = l; right = r; bottom = b; top = t;
  // Recompute distances
  if (wrapX) invDist.x = (right-left)/dX;
  else invDist.x = (right-left)/(dX-1);
  if (wrapY) invDist.y = (top-bottom)/dY;
  else invDist.y = (top-bottom)/(dY-1);
  LFactor = 2*(sqr(invDist.x)+sqr(invDist.y));
  invLFactor = 1./LFactor;
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
void FieldBase<T>::setWrap(bool x, bool y) {
  setWrapX(x); setWrapY(y);
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
void FieldBase<T>::setEdge(int edge, double value, bool lock) {
  switch (edge) {
  case 0: // Top
    for (int i=0; i<dX; i++) {
      at(i,dY-1) = value;
      lockAt(i,dY-1) = lock;
    }
    break;
  case 1: // Right
    for(int i=0; i<dY; i++) {
      at(dX-1,i) = value;
      lockAt(dX-1,i) = lock;
    }
    break;
  case 2: // Bottom
    for(int i=0; i<dX;i++) {
      at(i,0) = value;
      lockAt(i,0) = lock;
    }
    break;
  case 3: // Left
    for(int i=0; i<dY; i++) {
      at(0,i) = value;
      lockAt(0,i) = lock;
    }
    break;
  default: break; // Anything else
  }
}

template<typename T>
void FieldBase<T>::setAll(const T& value) {
  for (int i=0; i<dX*dY; i++) array[i] = value;
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
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds(x,y);
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
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds(x,y);
  return array[y*dX+x];
}

template<typename T>
T& FieldBase<T>::operator()(int x, int y) {
  return at(x,y);
}

template<typename T>
T FieldBase<T>::at(vect<> pos, bool thrw) const {
  correctPos(pos);
  if (!checkPos(pos, thrw)) return T(0);

  double X = (pos.x-left)*invDist.x, Y = (pos.y-bottom)*invDist.y;
  double bx = (int)X, by = (int)Y;

  double x = X-bx, y = Y-by;

  T tl = at(bx, by+1), tr = at(bx+1, by+1);
  T bl = at(bx, by), br = at(bx+1, by);
  T p_E = x*(br-bl)*invDist.x + bl;
  T p_F = x*(tr-tl)*invDist.x + tl;
  return y*(p_F-p_E)*invDist.y + p_E;
}

template<typename T>
T FieldBase<T>::operator()(vect<> pos, bool thrw) const {
  correctPos(pos);
  if (!checkPos(pos, thrw)) return T(0);

  double X = (pos.x-left)*invDist.x, Y = (pos.y-bottom)*invDist.y;
  double bx = (int)X, by = (int)Y;

  double x = X-bx, y = Y-by;

  T tl = at(bx, by+1), tr = at(bx+1, by+1);
  T bl = at(bx, by), br = at(bx+1, by);
  T p_E = x*(br-bl)*invDist.x +bl;
  T p_F= x*(tr-tl)*invDist.x +tl;
  return y*(p_F-p_E)*invDist.y + p_E;
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
string FieldBase<T>::printLocks() const {
  if (!usesLocks) return "0";
  stringstream stream;
  stream << '{';
  for (int y=dY-1; y>=0; y--) {
    stream << '{';
    for (int x=0; x<dX; x++) {
      stream << lockAt(x,y);
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
bool& FieldBase<T>::lockAt(int x, int y) {
  if (!usesLocks) throw NoLocks();
  if (wrapX) {
    while (x<0) x+=dX;
    while (x>=dX) x-=dX;
  }
  if (wrapY) {
    while (y<0) y+=dY;
    while (y>=dY) y-=dY;
  }
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds(x,y);
  
  return locks[x+dX*y];
}

template<typename T>
bool FieldBase<T>::lockAt(int x, int y) const {
  if (!usesLocks) throw NoLocks();
  if (wrapX) {
    while (x<0) x+=dX;
    while (x>=dX) x-=dX;
  }
  if (wrapY) {
    while (y<0) y+=dY;
    while (y>=dY) y-=dY;
  }
  if (x>=dX || x<0 || y>=dY || y<0) throw OutOfBounds(x,y);
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
    return 0.5*invDist.x*(at(x+1,y)-at(x-1,y));
    /*
    if (x==0) return 0.5*invDist.x*(at(1,y)-at(dX-1,y));
    else if (x==dX-1) return 0.5*invDist.x*(at(0,y)-at(dX-2,y));
    else return 0.5*invDist.x*(at(x+1,y)-at(x-1,y));
    */
  }
  else { // Don't wrap x
    if (x==0) return invDist.x*(at(1,y)-at(0,y));
    else if (x==dX-1) return invDist.x*(at(dX-1,y)-at(dX-2,y));
    else return 0.5*invDist.x*(at(x+1,y)-at(x-1,y));
  }
}

template<typename T>
T FieldBase<T>::DY(int x, int y) const {
  if (wrapY) {
    return 0.5*invDist.y*(at(x,y+1)-at(x,y-1));
    /*
    if (y==0) return 0.5*invDist.y*(at(x,1)-at(x,dY-1));
    else if (y==dY-1) return 0.5*invDist.y*(at(x,0)-at(x,dY-2));
    else return 0.5*invDist.y*(at(x,y+1)-at(x,y-1));
    */
  }
  else { // Don't wrap y
    if (y==0) return invDist.y*(at(x,1)-at(x,0));
    else if (y==dY-1) return invDist.y*(at(x,dY-1)-at(x,dY-2));
    else return 0.5*invDist.y*(at(x,y+1)-at(x,y-1));
  }
}

template<typename T>
T FieldBase<T>::D2X(int x, int y) const {
  if (wrapX) return sqr(invDist.x)*(at(x+1,y)-2*at(x,y)+at(x-1,y));
  else {
    if (x==0) return 2*D2X(1,y)-D2X(2,y); // Linearly interpolate
    else if (x==dX-1) return 2*D2X(dX-2,y)-D2X(dX-3,y); // Linearly interpolate
    else return sqr(invDist.x)*(at(x+1,y)-2*at(x,y)+at(x-1,y));
  }
}

template<typename T>
T FieldBase<T>::D2Y(int x, int y) const{
  if (wrapY) return sqr(invDist.y)*(at(x,y+1)-2*at(x,y)+at(x,y-1));
  else {
    if (y==0) return 2*D2X(x,1)-D2X(x,2); // Linearly interpolate
    else if (y==dY-1) return 2*D2X(x,dY-2)-D2X(x,dY-3); // Linearly interpolate
    else return sqr(invDist.y)*(at(x,y+1)-2*at(x,y)+at(x,y-1));
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
  int sx = wrapX?0:1, ex = wrapX?dX:dX-1;
  int sy = wrapY?0:1, ey = wrapY?dY:dY-1;
  for (int iter=0; iter<solveIterations && maxDelta>tollerance; iter++) {
    // Even squares
    maxDelta = 0;
    for (int y=sy; y<ey; y++)
      for (int x=sx; x<ex; x++)
        if ((x+y)%2==0 && (!usesLocks || !lockAt(x,y))) {
          T value = (1-omega)*at(x,y) + omega*invLFactor*(sqr(invDist.x)*(at(x+1,y)+at(x-1,y)) + sqr(invDist.y)*(at(x,y+1)+at(x,y-1)));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
    // Odd squares
    for (int y=sy; y<ey; y++)
      for (int x=sx; x<ex; x++)
	if ((x+y)%2 && (!usesLocks || !lockAt(x,y))) {
          T value = (1-omega)*at(x,y) + omega*invLFactor*(sqr(invDist.x)*(at(x+1,y)+at(x-1,y)) + sqr(invDist.y)*(at(x,y+1)+at(x,y-1)));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
    // Boundaries
    if (!wrapY)
      for (int x=0; x<dX; x++) {
	//if (!usesLocks || !lockAt(x,0)) at(x,0) = 2*at(x,1)-at(x,2);
        //if (!usesLocks || !lockAt(x,dY-1)) at(x,dY-1) = 2*at(x,dY-2)-at(x,dY-3);
        if (!usesLocks || !lockAt(x,0)) at(x,0) = at(x,1);
        if (!usesLocks || !lockAt(x,dY-1)) at(x,dY-1) = at(x,dY-2);

      }
    if (!wrapX)
      for (int y=0; y<dY; y++) {
        //if (!usesLocks || !lockAt(0,y)) at(0,y) = 2*at(1,y)-at(2,y);
	//if (!usesLocks || !lockAt(dX-1,y)) at(dX-1,y) = 2*at(dX-2,y)-at(dX-3,y);
	if (!usesLocks || !lockAt(0,y)) at(0,y) = at(1,y);
        if (!usesLocks || !lockAt(dX-1,y)) at(dX-1,y) = at(dX-2,y);
      }
  }
}

/// Successive Over-Relaxation (SOR) using checkerboard updating with source
template<typename T>
void FieldBase<T>::SOR_solver(FieldBase& source, double mult) {
  // Check that the source field has the same dimensions
  matches(&source);
  // Calculate relaxation parameter
  double omega = 2.0/(1+PI/dX);
  double maxDelta = 1e9;
  int sx = wrapX?0:1, ex = wrapX?dX:dX-1;
  int sy = wrapY?0:1, ey = wrapY?dY:dY-1;
  // Successively over-relax the system
  for (int iter=0; iter<solveIterations && maxDelta>tollerance; iter++) {
    maxDelta = 0;
    // Even squares
    for (int y=sy; y<ey; y++)
      for (int x=sx; x<ex; x++)
        if ((x+y)%2==0 && (!usesLocks || !lockAt(x,y))) {
          T value = (1-omega)*at(x,y) + omega*invLFactor*(sqr(invDist.x)*(at(x+1,y)+at(x-1,y)) + sqr(invDist.y)*(at(x,y+1)+at(x,y-1)) - mult*source(x,y));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
    // Odd squares
    for (int y=sy; y<ey; y++)
      for (int x=sx; x<ex; x++)
        if ((x+y)%2 && (!usesLocks || !lockAt(x,y))) {
          T value = (1-omega)*at(x,y) + omega*invLFactor*(sqr(invDist.x)*(at(x+1,y)+at(x-1,y)) + sqr(invDist.y)*(at(x,y+1)+at(x,y-1)) - mult*source(x,y));
          double delta = sqr(at(x,y)-value);
          if (delta>maxDelta) maxDelta = delta/sqr(at(x,y));
          at(x,y) = value;
        }
    // Boundaries
    if (!wrapY) 
      for (int x=0; x<dX; x++) {
	//if (!usesLocks || !lockAt(x,0)) at(x,0) = 2*at(x,1)-at(x,2);
	//if (!usesLocks || !lockAt(x,dY-1)) at(x,dY-1) = 2*at(x,dY-2)-at(x,dY-3);
	if (!usesLocks || !lockAt(x,0)) at(x,0) = at(x,1);
	if (!usesLocks || !lockAt(x,dY-1)) at(x,dY-1) = at(x,dY-2);

      }
    if (!wrapX) 
      for (int y=0; y<dY; y++) {
	//if (!usesLocks || !lockAt(0,y)) at(0,y) = 2*at(1,y)-at(2,y);
	//if (!usesLocks || !lockAt(dX-1,y)) at(dX-1,y) = 2*at(dX-2,y)-at(dX-3,y);
	if (!usesLocks || !lockAt(0,y)) at(0,y) = at(1,y);
	if (!usesLocks || !lockAt(dX-1,y)) at(dX-1,y) = at(dX-2,y);
      }
  }
}

template<typename T>
void FieldBase<T>::correctPos(vect<>& pos) const {
  double width = right-left, height = top-bottom;
  if (wrapX) {
    while (pos.x>right) pos.x -= width;
    while (pos.x<left) pos.x += width;
  }
  if (wrapY) {
    while (pos.y>top) pos.y -= height;
    while (pos.y<bottom) pos.y += height;
  }
}

template<typename T>
bool FieldBase<T>::checkPos(const vect<> pos, bool thrw) const {
  if (pos.x<left || right<pos.x || pos.y<bottom || top<pos.y) {
    if (thrw) throw InterpolateOutOfBounds(pos.x,pos.y);
    return false;
  }
  return true;
}

template<typename T>
template<typename S>
bool FieldBase<T>::matches(const FieldBase<S> *A) const{
  if (A->getDX()!=dX || A->getDY()!=dY) throw FieldMismatch();
}
