#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <sstream>
using std::stringstream;

double sqr(double x) { return x*x; }

/// Profile test class
class Profile {
public:
  Profile(int, double);
  void run(int);
  string print();
  
private:
  double vel(int);

  double rhoMax;
  double phi;
  double V;
  double epsilon;
  int size;
  double *array, *bb_array;
};

Profile::Profile(int width, double phi) {
  size = width;
  this->phi = phi;  
  array = new double[size];
  bb_array = new double[size];
  rhoMax = 1;
  epsilon = 0.01;
  V=1;

  for (int i=0; i<size; i++) array[i] = phi;
}

void Profile::run(int iters) {
  for (int j=0; j<size; j++) bb_array[j] = array[j];

  // Run the program
  for (int it=0; it<iters; it++) {
    // Update profile
    for (int j=1; j<size-1; j++) {
      double Vl = vel(j+1)-vel(j);
      double fluxL = (1-array[j-1]/rhoMax)*2*sqr(Vl)*array[j+1]*array[j]*epsilon;
      bb_array[j-1] += fluxL;
      double Vr = vel(j-1)-vel(j);
      double fluxR = (1-array[j+1]/rhoMax)*2*sqr(Vr)*array[j-1]*array[j]*epsilon;
      bb_array[j+1] += fluxR;
      bb_array[j] -= (fluxR+fluxL);
    }
    for (int j=0; j<size; j++) array[j] = bb_array[j];
  }
}

string Profile::print() {
  stringstream stream;
  stream << "{";
  double total = 0;
  for (int j=0; j<size; j++) total += array[j];
  for (int j=0; j<size; j++) {
    stream << "{" << j << "," << array[j] << "}";
    //stream << "{" << j << "," << array[j]/total << "}";
    if (j!=size-1) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

double Profile::vel(int i) {
  double R = (size-1)/2.;
  return V*(1-sqr((i-R)/R));
}


int main(int argc, char* argv[]) {

  double phi = 0.5;
  stringstream stream;
  if (argc>1) {
    stream << argv[1];
    stream >> phi;
  }
  cout << "phi=" << phi << ";" << endl;

  Profile profile(100, phi);
  
  profile.run(100000000); 
  cout << "L=" << profile.print() << ";\n";
  cout << "ListLinePlot[L,PlotStyle->Black,PlotRange->{0,1}]";
  
  return 0;
}
