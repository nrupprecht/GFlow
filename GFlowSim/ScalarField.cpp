#include "ScalarField.h"

ScalarField::ScalarField() : nsx(0), nsy(0), dx(0), idx(0), dy(0), idy(0), wrapX(false), wrapY(false), array(0), lap_array(0), diffusion(default_diffusion), lambda(default_lambda), printPoints(50) {};

ScalarField::ScalarField(double l, double r, double b, double t) : bounds(Bounds(l,r,b,t)), nsx(100), nsy(100), dx((r-l)/nsx), idx(1./dx), dy((t-b)/nsy), idy(1./dy), wrapX(false), wrapY(false), diffusion(default_diffusion), lambda(default_lambda), printPoints(50) {
    array = new double[nsx*nsy];
    lap_array = new double[nsx*nsy];
}

void ScalarField::set(std::function<double(double,double)> f) {
  double X, Y=bounds.bottom;
  for (int y=0; y<nsy; ++y) {
    X = bounds.left;
    for (int x=0; x<nsx; ++x) {
      array[nsx*y+x] = f(X,Y);
      X += dx;
    }
    Y += dy;
  }
}

void ScalarField::setBounds(Bounds b) {
  double res = (bounds.right-bounds.left)/nsx;
  bounds = b;
  if (0<res) setResolution(res);
}

void ScalarField::setBounds(double l, double r, double b, double t) {
  double res = (r-l)/nsx;
  bounds = Bounds(l,r,b,t);
  if (0<res) setResolution(res);
}

void ScalarField::setResolution(double res) {
  if (res<=0) throw IllegalResolution();
  if (array) delete [] array;
  if (lap_array) delete [] lap_array;

  nsx = (bounds.right - bounds.left)/res;
  nsy = (bounds.top - bounds.bottom)/res;
  dx = (bounds.right - bounds.left)/nsx;
  dy = (bounds.top - bounds.bottom)/nsy;
  ++nsx; ++nsy;
  idx = 1./dx; idy = 1./dy;
  array = new double[nsx*nsy];
  lap_array = new double[nsx*nsy];
  // Initialize values
  for (int i=0; i<nsx*nsy; ++i) array[i] = lap_array[i] = 0;
}

void ScalarField::set(Bounds b, double res) {
  bounds = b;
  setResolution(res);
}

void ScalarField::discard() {
  if (array) delete [] array;
  array = 0;
  if (lap_array) delete [] lap_array;
  lap_array = 0;
  nsx = nsy = 0; 
  dx = dy = idx = idy = 0;
}

double ScalarField::get(double x, double y) const {
  int X = (int)((x-bounds.left)*idx), Y = (int)((y-bounds.bottom)*idy);
  double cX = bounds.left + dx*X, cY = bounds.bottom +dy*Y;
  double DX = x-cX, DY = y-cY;
  double V1, V2;

  double A, B, C, D;

  if (-1<X && X+1<nsx && -1<Y && Y+1<nsy) {
    A = at(X,Y);
    B = at(X+1,Y);
    C = at(X,Y+1);
    D = at(X+1,Y+1);
  }
  else return 0;

  V1 = (B-A)*DX*idx + A;
  V2 = (D-C)*DX*idx + C;
  return (V2-V1)*DY*idy+V1;
}

double& ScalarField::at(int x, int y) {
  if (x<0 || nsx<=x) {
    if (wrapX && x<0) x+= nsx;
    else if (wrapX && nsx-1<x) x-=nsx;
    else throw ScalarField::FieldOutOfBounds(x,y);
  }
  if (y<0 || nsy<=y) {
    if (wrapY && y<0) y+=nsy;
    else if (wrapY && nsy-1<y) y-=nsy;
    else throw ScalarField::FieldOutOfBounds(x,y);
  }
  return array[nsx*y+x];
}

double ScalarField::at(int x, int y) const {
  if (x<0 || nsx<=x) {
    if (wrapX && x<0) x+= nsx;
    else if (wrapX && nsx-1<x) x-=nsx;
    else throw ScalarField::FieldOutOfBounds(x,y);
  }
  if (y<0 || nsy<=y) {
    if (wrapY && y<0) y+=nsy;
    else if (wrapY && nsy-1<y) y-=nsy;
    else throw ScalarField::FieldOutOfBounds(x,y);
  }
  return array[nsx*y+x];
}

double ScalarField::lap(int x, int y) const {
  if (x<0 || nsx<=x || y<0 || nsy<=y) throw ScalarField::FieldOutOfBounds(x,y);
  return lap_array[nsx*y+x];
}

double ScalarField::dX(double x, double y) const {
  int X = (x-bounds.left)*idx, Y = (y-bounds.bottom)*idy;
  return dXAt(X,Y);
}

double ScalarField::dY(double x, double y) const {
  int X = (x-bounds.left)*idx, Y = (y-bounds.bottom)*idy;
  return dYAt(X,Y);
}

double ScalarField::dXAt(int x, int y) const {
  if (wrapX || (0<x && x<nsx-1)) return 0.5*(at(x+1, y)-at(x-1,y))*idx;
  else if (x==0)      return (at(1,y)-at(0,y))*idx;
  else /* x==nsx-1 */ return 0.5*(at(nsx-1, y)-at(nsx-2, y))*idx;
}

double ScalarField::dYAt(int x, int y) const {
  if (wrapY || (0<y && y<nsy-1)) return 0.5*(at(x,y+1)-at(x,y-1))*idy;
  else if (y==0)      return(at(x,1)-at(x,0))*idy;
  else /* y==nsy-1 */ return 0.5*(at(x,nsy-1)-at(x,nsy-2))*idy;
}

double ScalarField::d2X(double x, double y) const {
  int X = (x-bounds.left)*idx, Y = (y-bounds.bottom)*idy;
  return d2XAt(X,Y);
}

double ScalarField::d2Y(double x, double y) const {
  int X = (x-bounds.left)*idx, Y = (y-bounds.bottom)*idy;
  return d2YAt(X,Y);
}

double ScalarField::d2XAt(int x, int y) const {
  if (wrapX || (0<x && x<nsx-1)) return (at(x+1,y)+at(x-1,y)-2*at(x,y))*sqr(idy);
  else if (x==0)      return 2*d2XAt(1,y)-d2XAt(2,y); 
  else /* x==nsx-1 */ return 2*d2XAt(nsx-2,y)-d2XAt(nsx-3,y);
}

double ScalarField::d2YAt(int x, int y) const {
  if (wrapY || (0<y && y<nsy-1)) return (at(x,y+1)+at(x,y-1)-2*at(x,y))*sqr(idy);
  else if (y==0)      return 2*d2YAt(x,1)-d2YAt(x,2);
  else /* y==nsy-1 */ return 2*d2YAt(x,nsy-2)-d2YAt(x,nsy-3);
}

vector<double> ScalarField::sliceX(double y, int res) {
  if (y<bounds.bottom || bounds.top<y) return vector<double>();
  vector<double> array;
  double DX = (bounds.right-bounds.left)/res;
  double X = bounds.left;
  for (int i=0; i<res; ++i) {
    array.push_back(get(X,y));
    X += DX;
  }
  return array;
}

vector<double> ScalarField::sliceY(double x, int res) {
  if (x<bounds.left || bounds.right < x) return vector<double>();
  vector<double> array;
  double DY = (bounds.top-bounds.bottom)/res;
  double Y = bounds.bottom;
  for (int i=0; i<res; ++i) {
    array.push_back(get(x,Y));
    Y += DY;
  }
  return array;
}

/*
double ScalarField::derivative(double x, double y, int nx, int ny) const {
  int X = (x-bounds.left)*idx, Y = (y-bounds.bottom)*idy;
  return derivativeAt(X, Y, nx, ny);
}

double ScalarField::derivativeAt(int x, int y, int nx, int ny) const {
  if (nx==0) {
    if (ny==0) return at(x,y);
    else { // Take y derivative
      
      if (0<y && y<nsy-1) return 0.5*(derivativeAt(x,y+1,nx,ny-1)-derivativeAt(x,y-1,nx,ny-1))*idy;
      else if (y==0)      return (derivativeAt(x,1,nx,ny-1)-derivativeAt(x,0,nx,ny-1))*idy;
      else                return 0.5*(derivativeAt(x,nsy-1,nx,ny-1)-derivativeAt(x,nsy-2,nx,ny-1))*idy;
    }
  }
  else { // nx != 0
    if (0<x && x<nsx-1) return 0.5*(derivativeAt(x+1,y,nx-1,ny)-derivativeAt(x-1,y,nx-1,ny))*idx;
    else if (x==0)      return (derivativeAt(1,y,nx-1,ny)-derivativeAt(0,y,nx-1,ny))*idx;
    else                return 0.5*(derivativeAt(nsx-1,y,nx-1,ny)-derivativeAt(nsx-2,y,nx-1,ny))*idx;
  }
}
*/

void ScalarField::laplacian() {
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x)
      lap_array[nsx*y+x] = d2XAt(x,y) + d2YAt(x,y);
}

void ScalarField::laplacian(ScalarField& sf) {
  sf.discard();
  sf.nsx = nsx; sf.nsy = nsy;
  sf.dx = dx; sf.idx = idx; sf.dy = dy; sf.idy = idy;
  if (nsx*nsy==0) sf.array = sf.lap_array = 0;
  else {
    sf.array = new double[nsx*nsy];
    sf.lap_array = new double[nsx*nsy];
  }
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x)
      sf.array[nsx*y+x] = d2XAt(x,y) + d2YAt(x,y);
}

void ScalarField::increase(double x, double y, double c) {
  if (!bounds.contains(x,y)) return;
  // Find field coordinates
  int X = (x-bounds.left)*idx, Y = (y-bounds.bottom)*idy;
  
  cinc(X,Y,c);
  return;

  // Two point gaussian
  cinc(X-2,Y+2,0.13*c); cinc(X-1,Y+2,0.22*c); cinc(X,Y+2,0.37*c); cinc(X+1,Y+2,0.22*c); cinc(X+2,Y+2,0.13*c);
  cinc(X-2,Y+1,0.23*c); cinc(X-1,Y+1,0.60*c); cinc(X,Y+1,0.78*c); cinc(X+1,Y+1,0.60*c); cinc(X+2,Y+1,0.23*c);
  cinc(X-2,Y,0.37*c);   cinc(X-1,Y,0.78*c);   cinc(X,Y,c);        cinc(X+1,Y,0.78*c);   cinc(X+2,Y,0.37*c);
  cinc(X-2,Y-1,0.23*c); cinc(X-1,Y-1,0.60*c); cinc(X,Y-1,0.78*c); cinc(X+1,Y-1,0.60*c); cinc(X+2,Y-1,0.23*c);
  cinc(X-2,Y-2,0.13*c); cinc(X-1,Y-2,0.22*c); cinc(X,Y-2,0.37*c); cinc(X+1,Y-2,0.22*c);cinc(X+2,Y-2,0.13*c);
}

void ScalarField::reduce(double x, double y, double change) {
  if (!bounds.contains(x,y)) return;
  // Find field coordinates
  int X = (x-bounds.left)*idx, Y = (y-bounds.bottom)*idy;
  at(X,Y) += change;
}

void ScalarField::propReduceCue(double x, double y, double lam) {
  if (!bounds.contains(x,y)) return;
  // Find field coordinates
  int X = (x-bounds.left)*idx, Y = (y-bounds.bottom)*idy;
  cue(X,Y) += lam;
}


void ScalarField::propReduceExec() {
  for (const auto &p : propReduce) {
    int x = p.first%nsx, y = p.first/nsx;
    at(x,y) += (p.second*at(x,y));
  }
  propReduce.clear();
}

void ScalarField::update(double epsilon) {  
  if (array==0) return;
  // Take the laplacian
  if (diffusion!=0) { // Diffusion
    laplacian();
    // Update field
    if (lambda!=0) // Decay
      for (int y=0;y<nsy; ++y)
	for(int x=0; x<nsx; ++x)
	  at(x,y) += (diffusion*lap(x,y) + lambda*at(x,y))*epsilon;
    else // No decay
      for (int y=0;y<nsy; ++y)
	for(int x=0; x<nsx; ++x)
	  at(x,y) += diffusion*lap(x,y)*epsilon;
  }
  else { // No diffusion
    if (lambda!=0) // Decay
      for (int y=0;y<nsy; ++y)
	for(int x=0; x<nsx; ++x)
	  at(x,y) += lambda*at(x,y)*epsilon;
    else; // Nothing needs to be done to update the field
  }
}

std::ostream& operator<<(std::ostream& out, const ScalarField& field) {
  if (field.array==0) {
    out << "{}";
    return out;
  }
  Bounds bounds = field.bounds;
  int printPoints = field.printPoints;
  out << "{";
  double X, Y = bounds.bottom, DX, DY;
  // Set step size
  if (printPoints>0) {
    DX = (bounds.right-bounds.left)/printPoints;
    DY = (bounds.top-bounds.bottom)/printPoints;
  }
  else {
    DX = field.dx;
    DY = field.dy;
  }
  if (printPoints>0) {
    for (int y=0; y<printPoints; ++y) {
      X = bounds.left;
      for (int x=0; x<printPoints; ++x) {
        out << "{" << X << "," << Y << "," << field.get(X,Y) << "}";
        X += DX;
	if (printPoints*y+x!=sqr(printPoints)-1) out << ',';
      }
      Y += DY;
    }
  }
  else { // Use full resolution
    int nsx = field.nsx, nsy = field.nsy;
    for (int y=0; y<nsy; ++y) {
      X = bounds.left;
      for (int x=0; x<nsx; ++x) {
        out << "{" << X << "," << Y << "," << field.array[nsx*y+x] << "}";
        X += DX;
	if (nsx*y+x!=nsx*nsy-1) out << ',';
      }
      Y += DY;
    }
  }
  out << "}";
  return out;
}

bool ScalarField::printToCSV(string filename) const {
  if (array==0) return true;
  if (max(nsx,nsy)<printPoints) printPoints=max(nsx,nsy);
  // Open file
  ofstream fout(filename);
  if(fout.fail()) return false;
  // Print data
  double X, Y = bounds.bottom, DX, DY;
  // Set step size
  if (printPoints>0) {
    DX = (bounds.right-bounds.left)/printPoints;
    DY = (bounds.top-bounds.bottom)/printPoints;
  }
  else { 
    DX = dx; 
    DY = dy;
  }
  if (printPoints>0) 
    for (int y=0; y<printPoints; ++y) {
      X = bounds.left;
      for (int x=0; x<printPoints; ++x) {
	fout << X << "," << Y << "," << get(X,Y) << endl;
	X += DX;
      }
      Y += DY;
    }
  else // Use full resolution
    for (int y=0; y<nsy; ++y) {
      X = bounds.left;
      for (int x=0; x<nsx; ++x) {
        fout << X << "," << Y << "," << array[nsx*y+x] << endl;
	X += DX;
      }
      Y += DY;
    }
  return true;
}

bool ScalarField::loadFromCSV(string filename) {
  ifstream fin(filename);
  if (fin.fail()) return false;
  // Read in values - they should be in the form [x,y,v\n]
  vector<Trio> values;
  double x, y, v;
  char c; // For the ','
  while (!fin.eof()) {
    fin >> x >> c >> y >> c >> v;
    values.push_back(Trio(x,y,v));
  }
  if (values.empty()) {
    discard(); // Empty field
    return true; 
  }
  // Find bounds and resolution
  double minX=1e9, maxX=-1e9, minY=1e9, maxY=1e-9;
  double DX=1e9, DY=1e9;
  Trio V = *values.begin();
  for (auto t : values) {
    // Look for bounds
    if (t.x<minX) minX = t.x;
    if (maxX<t.x) maxX = t.x;
    if (t.y<minY) minY = t.y;
    if (maxY<t.y) maxY = t.y;
    // Look for resolution
    double tdx = fabs(V.x-t.x), tdy = fabs(V.y-t.y);
    if (tdx<DX && 0<tdx) DX = fabs(V.x-t.x);
    if (tdy<DY && 0<tdy) DY = fabs(V.y-t.y);
  }
  
  // Set values
  bounds.left = minX; bounds.right = maxX;
  bounds.bottom = minY; bounds.top = maxY;
  dx = DX; dy = DY;
  nsx = (bounds.right - bounds.left)/dx+1;
  nsy = (bounds.top - bounds.bottom)/dy+1;
  idx = 1./dx; idy = 1./dy;
  array = new double[nsx*nsy];
  lap_array = new double[nsx*nsy];
  // Set values
  for (const auto& v : values) {
    int x = (v.x-bounds.left)*idx;
    int y = (v.y-bounds.bottom)*idy;
    array[nsx*y+x] = v.z;
  }
  // Return success
  return true;
}

void ScalarField::cinc(int x, int y, double increase) {
  if (-1<x && x<nsx && -1<y && y<nsy) at(x,y) += increase;
}

double& ScalarField::cue(int x, int y) {
  auto p = propReduce.find(nsx*y+x);
  if (p!=propReduce.end()) return p->second;
  else {
    propReduce.insert(pair<int,double>(nsx*y+x,0));
    return propReduce.find(nsx*y+x)->second;
  }
}
