#include "Fluid.h"

Fluid::Fluid() : dX(10), dY(10), epsilon(0.001), time(0), buffer(0), bbuffer(1), wrap(false) {
  gravity = vect<>(0,-1);
  // Set up array
  array[0] = new Element*[dY];
  array[1] = new Element*[dY];
  for (int i=0; i<dY; i++) {
    array[0][i] = new Element[dX];
    array[1][i] = new Element[dX];
  }
  /*
  for (int x=0; x<dX; x++) {
    array[buffer][0][x].rho = 0.1;
    array[buffer][dY-1][x].rho = 0.1;
  }
  for (int y=0; y<dY; y++) {
    array[buffer][y][0].rho = 0.1;
    array[buffer][y][dX-1].rho = 0.1;
  }
  */
}

Fluid::~Fluid() {
  for (int i=0; i<dX; i++) delete [] array[0][i];
  delete [] array[0];

  for (int i=0; i<dX; i++) delete [] array[1][i];
  delete [] array[1];
}

void Fluid::run(double runTime) {
  time = 0;
  while (time<runTime) {
    // Calculate body forces
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
	forces(x,y);
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
	conditionF(x,y);

    //------>
    if (true) {
      cout << "Forces:\n";
      for (int y=dY-1; y>=0; y--) {
	for (int x=0; x<dX; x++)
	  cout << array[buffer][y][x].force << '\t';
	cout << endl;
      }
      cout << endl;
    }
    // ------<

    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
	pressures(x,y);
    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
	conditionP(x,y);

    //------>
    if (true) {
      cout << "Pressure:\n";
      for (int y=dY-1; y>=0; y--) {
	for (int x=0; x<dX; x++)
	  cout << array[bbuffer][y][x].pressure << '\t';
	cout << endl;
      }
      cout << endl;
    }
    // ------<

    for (int y=0; y<dY; y++)
      for (int x=0; x<dX; x++)
        updates(x,y);

    //------>
    if (true) {
      cout << "Rho:\n";
      for (int y=dY-1; y>=0; y--) {
	for (int x=0; x<dX; x++)
	  cout << array[bbuffer][y][x].rho << '\t';
	cout << endl;
      }
      cout << endl;
    }
    // ------<

    //------>
    if (false) {
      cout << "Velocity:\n";
      for (int y=dY-1; y>=0; y--) {
	for (int x=0; x<dX; x++)
	  cout << array[bbuffer][y][x].fU << '\t';
	cout << endl;
      }
      cout << endl;
    }
    // ------<

    // Swap buffers
    swap(buffer, bbuffer);

    // Update time
    time += epsilon;
  }
}

// Find the net forces on a fluid element
inline void Fluid::forces(int x, int y) {
  vect<> force;
  // Gravity
  force = rhoAt(x,y)*gravity;
  // Pressure from surrounding elements
  if (y-1>=0) force += array[buffer][y-1][x].pressure*E1;
  if (y+1<dY) force += array[buffer][y+1][x].pressure*(-1*E1);
  if (x-1>=0) force += array[buffer][y][x-1].pressure*E0;
  if (x+1<dX) force += array[buffer][y][x+1].pressure*(-1*E0);
  // ---> Calculate fluid - solid forces here

  // Boundary conditions
  if (y==0) force = pclamp(force, E1);
  else if (y==dY-1) force = pclamp(force, -E1);
  else if (x==0) force = pclamp(force, E0);
  else if (x==dX-1) force = pclamp(force, -E0);

  forceAt(x,y) = force;
}

inline void Fluid::conditionF(int x, int y) {
  if (x==0) forceAt(0,y) = -pclamp(forceAt(1,y), -E0);
  else if (x==dX-1) forceAt(dX-1,y) = -pclamp(forceAt(dX-2,y), E0);
  else if (y==0) forceAt(x,0) =-pclamp(forceAt(x,1), -E1);
  else if (y==dY-1) forceAt(x,dY-1) = -pclamp(forceAt(x,dY-2), E1);  
}

// Calculate the pressure on a fluid element
inline void Fluid::pressures(int x, int y) {
  // Reset
  bbPAt(x,y) = 0;
  // Calculate pressure from inward normal forces
  int count = 0;
  if (y>0) bbPAt(x,y) += clamp(forceAt(x,y-1)*vect<>(0,1)) + pAt(x,y-1);
  if (y+1<dY) bbPAt(x,y) += clamp(forceAt(x,y+1)*vect<>(0,-1)) + pAt(x,y+1);
  if (x>0) bbPAt(x,y) += clamp(forceAt(x-1,y)*vect<>(1,0)) + pAt(x-1,y);
  if (x+1<dX) bbPAt(x,y) += clamp(forceAt(x+1,y)*vect<>(-1,0)) + pAt(x+1,y);
  bbPAt(x,y) *= 0.25;
}

inline void Fluid::conditionP(int x, int y) {
  if (x==0) {
    pAt(x,y) = pAt(x+1,y);
    bbPAt(x,y) = pAt(x+1,y);
  }
  else if (x==dX-1) {
    pAt(x,y) = pAt(x-1,y);
    bbPAt(x,y) = pAt(x-1,y);
  }
  else if (y==dY-1) {
    pAt(x,y) = 0;
    bbPAt(x,y) = 0;
  }
  else if (y==0 || y==dY-1) {
    pAt(x,y) = pAt(x,y+1);
    bbPAt(x,y) = pAt(x,y+1);
 }
}

inline void Fluid::updates(int x, int y) {
  vect<> delrho = delRho(x,y);
  vect<> advection = advectU(x,y);
  vect<> delp = delP(x,y);
  double divu = divU(x,y);
  // Calculate new rho
  //double delta = array[bbuffer][y][x].fU*delrho + array[bbuffer][y][x].rho*divu;
  //array[bbuffer][y][x].rho = array[buffer][y][x].rho - epsilon*delta;

  // Calculate new U
  vect<> dV = advection + (1/array[buffer][y][x].rho)*delp - /*nu*delSqrU -*/ array[buffer][y][x].force; //**
  vect<> velocity = array[buffer][y][x].fU - epsilon*dV;

  // Boundary conditions
  if (y==0) velocity = pclamp(velocity, E1);
  else if (y==dY-1) velocity = pclamp(velocity, -E1);
  else if (x==0) velocity = pclamp(velocity, E0);
  else if (x==dX-1) velocity = pclamp(velocity, -E0);
  array[bbuffer][y][x].fU = velocity;
}

vect<> Fluid::dU_dx(int x, int y) {
  if (x+1<dX && x>0) return 0.5*(array[buffer][y][x+1].fU - array[buffer][y][x-1].fU);
  else if (x==0) return (array[buffer][y][x+1].fU - array[buffer][y][x].fU);
  else return (array[buffer][y][x].fU - array[buffer][y][x-1].fU);
}

vect<> Fluid::dU_dy(int x, int y) {
  if (y+1<dY && y>0) return 0.5*(array[buffer][y+1][x].fU - array[buffer][y-1][x].fU);
  else if (y==0) return (array[buffer][y+1][x].fU - array[buffer][y][x].fU);
  else return (array[buffer][y][x].fU - array[buffer][y-1][x].fU);
}

vect<> Fluid::delRho(int x, int y) {
  double dRho_dx, dRho_dy;
  if (x+1<dX && x>0) dRho_dx = 0.5*(array[buffer][y][x+1].rho - array[buffer][y][x-1].rho);
  else if (wrap) { // wrapping
    if (x==0) dRho_dx = (rhoAt(x+1,y)- rhoAt(dX-1,y));
    else dRho_dx = (rhoAt(0,y) - rhoAt(dX-2,y));
  }
  else { // No wrapping
    if (x==0) dRho_dx = (rhoAt(x+1,y)- rhoAt(x,y));
    else dRho_dx = (rhoAt(x,y) - rhoAt(x-1,y));
  }
  if (y+1<dY && y>0) dRho_dx = 0.5*(rhoAt(x,y+1) - rhoAt(x,y-1));
  else if (wrap) {
    if (y==0) dRho_dy = (rhoAt(x,1) - rhoAt(x,dY-1));
    else dRho_dy = (rhoAt(x,0) - rhoAt(x,y-1));
  }
  else { // No wrapping
    if (y==0) dRho_dy = (rhoAt(x,y+1) - rhoAt(x,y));
    else dRho_dy = (rhoAt(x,y) - rhoAt(x,y-1));
  }
  return vect<>(dRho_dx, dRho_dy);
}

vect<> Fluid::delP(int x, int y) {
  double dP_dx, dP_dy;
  if (x+1<dX && x>0) dP_dx = 0.5*(array[buffer][y][x+1].pressure - array[buffer][y][x-1].pressure);
  else if (x==0) dP_dx = (array[buffer][y][x+1].pressure - array[buffer][y][x].pressure);
  else dP_dx = (array[buffer][y][x].pressure - array[buffer][y][x-1].pressure);
  if (y+1<dY && y>0) dP_dx = 0.5*(array[buffer][y+1][x].pressure - array[buffer][y-1][x].pressure);
  else if (y==0) dP_dy = (array[buffer][y+1][x].pressure - array[buffer][y][x].pressure);
  else dP_dy = (array[buffer][y][x].pressure - array[buffer][y-1][x].pressure);

  return vect<>(dP_dx, dP_dy);
}

double Fluid::divU(int x, int y) {
  return dU_dx(x,y).x + dU_dy(x,y).y;
}

vect<> Fluid::advectU(int x, int y) {
  vect<> V = array[buffer][y][x].fU;
  vect<> dudx = dU_dx(x,y), dudy = dU_dy(x,y);
  return vect<>(V.x*dudx.x+V.y*dudy.x, V.x*dudx.y+V.y*dudy.y);
}
