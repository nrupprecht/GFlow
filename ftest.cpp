#include "GFlow.h"

int main(int argc, char* argv[]) {
  int dim = 25;
  double time = 0.1;

  // See rands
  srand48(std::time(0));
  srand(std::time(0));

  stringstream stream;
  if (argc>1) {
    stream << argv[1];
    stream >> dim;
  }

  if (argc>2) {
    stream.clear();
    stream << argv[2];
    stream >> time;
  }

  GFlow test(dim,dim);

  //***** Test if this works
  // test.setResolution(2*dim,2*dim);
  // test.setBounds(0,2,0,2);
  //*****

  test.setGravity(vect<>());
  // test.setViscosity(0.1);
  // test.setEpsilon(1e-4);

  /*
  double R = 0.1;
  Particle *P = new Particle(vect<>(0.5, 0.85), R);
  P->setMass(0.1);
  test.addParticle(P);

  P = new Particle(vect<>(0.2, 0.1), R);
  test.addParticle(P);

  P = new Particle(vect<>(0.8, 0.1), R);
  test.addParticle(P);

  test.addWall(new Wall(vect<>(0,0), vect<>(1,0), true));
  test.addWall(new Wall(vect<>(0,1), vect<>(1,1), true));
  */

  //test.addWall(new Wall(vect<>(0,0), vect<>(1,0)));

  test.setDispDelay(1./60.);
  test.run(time);
  cout << "Epsilon: " << test.getEpsilon() << endl;
  cout << "Time: " << test.getRealTime() << endl;
  
  /*
  cout << "Pbdd=" << test.printP_bdd() << ";\n";
  cout << "Ubdd=" << test.printU_bdd() << ";\n";
  cout << "Vbdd=" << test.printV_bdd() << ";\n";
  cout << "coeff=" << test.printC() << ";\n";
  */

  /*
  cout << "press=" << test.getPressureRec() << ";\n";
  cout << test.printPositionRec() << "\n";
  cout << "walls=" << test.printWalls() << ";\n";
  cout << "R=" << test.printRadiusRec() << ";\n";
  cout << test.printPositionAnimationCommand("frames1") << endl;
  cout << test.printPressureAnimationCommand(true, "press","frames2") << endl;
  cout << "frames3=Table[Show[frames2[[i]],frames1[[i]]],{i,1,Length[frames1]}];\n";
  cout << "ListAnimate[frames3]";
  */

  cout << "p=" << test.printPressure() << ";\n";
  cout << "p3d=" << test.printPressure3D() << ";\n";
  cout << "vfmag=" << test.printVFNorm() << ";\n";
  cout << "vf=" << test.printVF() << ";\n";
  cout << "vfn=" << test.printVFN() << ";";
  
  return 0;

}

// For dt = 0.001, nx=ny=50 is the max
