#include "GFlow.h"

int main(int argc, char* argv[]) {
  int dim = 80;
  double time = 0.001;

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

  GFlow test(2*dim,dim);
  test.setBounds(0,2,0,1);
  test.setViscosity(1);

  test.addWall(new Wall(vect<>(0.5,0.75), vect<>(1.5,0.75)));
  test.addWall(new Wall(vect<>(0.5,0.5),vect<>(1.5,0.5)));
  test.addWall(new Wall(vect<>(0.5,0.25),vect<>(1.5,0.25)));
  test.addParticle(new Particle(vect<>(0.2,0.2),0.1));
  test.setGravity(Zero);

  for (int i=1; i<=dim; i++) {
    test.lockP(1,i,100);
    test.lockP(2*dim,i,-100);
  }

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

  cout << "press=" << test.getPressureRec() << ";\n";
  cout << test.printPositionRec() << "\n";
  cout << "walls=" << test.printWalls() << ";\n";
  cout << "R=" << test.printRadiusRec() << ";\n";
  cout << test.printPositionAnimationCommand("frames1") << endl;
  cout << test.printPressureAnimationCommand(true, "press","frames2") << endl;
  cout << "frames3=Table[Show[frames2[[i]],frames1[[i]]],{i,1,Length[frames1]}];\n";
  cout << "ListAnimate[frames3]" << endl;
  cout << "vel=" << test.getVelocityRec() << ";\n";
  cout << test.printVFAnimationCommand() << endl;


  cout << "p=" << test.printPressure() << ";\n";
  cout << "p3d=" << test.printPressure3D() << ";\n";
  cout << "vfmag=" << test.printVFNorm() << ";\n";
  cout << "vf=" << test.printVF() << ";\n";
  cout << "vfn=" << test.printVFN() << ";";

  return 0;

}
