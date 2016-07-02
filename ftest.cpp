#include "GFlow.h"

int main(int argc, char* argv[]) {
  int dim = 50;
  double time = 1;

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
  //test.addParticles(1, 0.1, 0, 0, 1, 0, 1);
  test.addParticle(new Particle(vect<>(0.5, 0.7), 0.2));
  test.addWall(new Wall(vect<>(0,0), vect<>(1,0), true));
  test.setDispDelay(0.005);
  test.run(time);

  cout << "Time: " << test.getRealTime() << endl;
  /*
  cout << "Iters: " << iter << endl;
  cout << "p=" << printPressure() << ";\n";
  cout << "p3d=" << printPressure3D() << ";\n";
  cout << "vfn=" << printVFNorm() << ";\n";
  cout << "vf=" << printVF() << ";\n";
  */

  cout << "Press=" << test.getPressureRec() << ";";

  return 0;

}

// For dt = 0.001, nx=ny=50 is the max
