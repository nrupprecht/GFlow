#include "ScalarField.h"

ScalarField::ScalarField() : nsx(0), nsy(0), dx(0), idx(0), dy(0), idy(0), array(0), lap_array(0) {};

ScalarField::ScalarField(double l, double r, double b, double t) : bounds(Bounds(l,r,b,t)), nsx(100), nsy(100), dx((r-l)/nsx), idx(1./dx), dy((t-b)/nsy), idy(1./dy) {
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

double& ScalarField::at(int x, int y) {
  if (x<0 || nsx<=x || y<0 || nsy<=y) throw ScalarField::FieldOutOfBounds(x,y);
  return array[nsx*y+x];
}

double ScalarField::at(int x, int y) const {
  if (x<0 || nsx<=x || y<0 || nsy<=y) throw ScalarField::FieldOutOfBounds(x,y);
  return array[nsx*y+x];
}

double ScalarField::lap(int x, int y) const {
  if (x<0 || nsx<=x || y<0 || nsy<=y) throw ScalarField::FieldOutOfBounds(x,y);
  return lap_array[nsx*y+x];
}

double ScalarField::dX(int x, int y) const {
  if (0<x && x<nsx-1) return 0.5*(at(x+1, y)+at(x-1, y)-2*at(x,y))*idx;
  else if (x==0)      return (at(1,y)-at(0,y))*idx;
  else /* x==nsx-1 */ return 0.5*(at(nsx-1, y)-at(nsx-2, y))*idx;
}

double ScalarField::dY(int x, int y) const {
  if (0<y && y<nsy-1) return 0.5*(at(x,y+1)+at(x,y-1)-2*at(x,y))*idy;
  else if (y==0)      return(at(x,1)-at(x,0))*idy;
  else /* y==nsy-1 */ return 0.5*(at(x,nsy-1)-at(x,nsy-2))*idy;
}

double ScalarField::derivative(int x, int y, int nx, int ny) const {
  if (nx==0) {
    if (ny==0) return at(x,y);
    else { // Take y derivative
      
      if (0<y && y<nsy-1) return 0.5*(derivative(x,y+1,nx,ny-1)+derivative(x,y-1,nx,ny-1)-2*derivative(x,y,nx,ny-1))*idy;
      else if (y==0)      return (derivative(x,1,nx,ny-1)-derivative(x,0,nx,ny-1))*idy;
      else /* y==nsy-1 */ return 0.5*(derivative(x,nsy-1,nx,ny-1)-derivative(x,nsy-2,nx,ny-1))*idy;
    }
  }
  else { // nx != 0
    if (0<x && x<nsx-1) return 0.5*(derivative(x+1,y,nx-1,ny)+derivative(x-1,y,nx-1,ny)-2*derivative(x,y,nx-1,ny))*idx;
    else if (x==0)      return (derivative(1,y,nx-1,ny)-derivative(0,y,nx-1,ny))*idx;
    else /* x==nsx-1 */ return 0.5*(derivative(nsx-1,y,nx-1,ny)-derivative(nsx-2,y,nx-1,ny))*idx;
  }
  
}

void ScalarField::laplacian() {
  for (int y=0; y<nsy; ++y)
    for (int x=0; x<nsx; ++x)
      lap_array[nsx*y+x] = derivative(x,y,2,0) + derivative(x,y,0,2);
}

void ScalarField::update(double epsilon, double diffusion, double lambda) {
  if (array==0) return;
  // Take the laplacian
  laplacian();
  // Update field
  for (int y=0;y<nsy; ++y)
    for(int x=0; x<nsx; ++x)
      at(x,y) += (diffusion*lap(x,y) - lambda*at(x,y))*epsilon;
}

std::ostream& operator<<(std::ostream& out, ScalarField& field) {
  out << "{";
  int nsx = field.nsx, nsy = field.nsy;
  double X, Y = field.bounds.bottom;
  for (int y=0; y<nsy; ++y) {
    X = field.bounds.left;
    for (int x=0; x<nsx; ++x) {
      out << "{" << X << "," << Y << "," << field.array[nsx*y+x] << "}";
      if (x!=nsx-1) out << ",";
      X += field.dx;
    }
    if (y!=nsy-1) out << ",";
    Y += field.dy;
  }
  out << "}";
  return out;
}
