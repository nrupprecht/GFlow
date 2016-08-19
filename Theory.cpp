#include "Theory.h"

Theory::Theory(int size, double sigma, double phi) : array(0), d_array(0) {
  propIters = 20;
  epsilon = 0.001;
  drag = 0;
  V0 = 0.5;
  R0 = 0.5;
  runTime = 0;
  
  initialize(size, sigma, phi);
  initialCondition();
}

void Theory::initialize(int size, double sigma, double phi) {
  this->size = size;
  this->sigma = sigma;
  this->phi = phi;

  scale = 2*R0/(double)size;
  cutoff = 2*sigma/scale;

  if (array) delete [] array;
  if (d_array) delete [] d_array;

  array = new double[size];
  d_array = new double[size];
  for (int i=0; i<size; i++) array[i] = d_array[i] = 0;
}

void Theory::solve(int max) {

  cout << "Radius: " << sigma << "\n";
  cout << "Cutoff: " << cutoff << " bins;\n";
  cout << "Scale: " << scale << ";\n";
  cout << "Initial: " << phi/(PI*sqr(sigma)) << ";\n";

  int i;
  double *deposit = new double[size];

  /// Do at most [max] iters
  auto start_t = clock();
  for (int it=0; it<max; it++) {
    // Iterate through bins
    for (int i=0; i<size; i++) d_array[i] = 0;
    for (int rr=1; rr<size-1; rr++) {
      double Fl = freqL(rr), Fr = freqR(rr);
      double frequency = Fl+Fr;
      Fl *= epsilon; Fr *= epsilon;
      
      /// Objects that are knocked out of the bin
      d_array[rr] -= frequency*epsilon;

      /// Knocked objects settle in other bins
      getPropagation(rr, Fl, Fr, deposit);
      for (int i=0; i<size; i++) d_array[i] += deposit[i];
    }
    // Every bin has been calculated, update main buffer
    for (int i=0; i<size; i++) array[i] += d_array[i];

    // Monitor change
    double diff = int_d_array();
    //if (diff>1e9) exit(0);
    //if (it%10==0 || diff==0) cout << "Diff: " << diff << endl; 

    if (diff<1e-4) break; // End condition
  }  
  auto end_t = clock();
  runTime = (double)(end_t-start_t)/CLOCKS_PER_SEC;
  delete [] deposit;
}

void Theory::test(string dist, double sigma, int iters, double V0) {

  this->sigma = sigma;
  this->V0 = V0;

  // Extract dist
  vector<double> vals;
  string str;
  double v;
  int i=0;
  for (auto c : dist) {
    if (c=='{' || c=='}');
    else if (c==',') {
      stringstream stream;
      stream << str;
      str.clear();
      stream >> v;
      vals.push_back(v);
    }
    else str += c;
  }
  // Get the last number
  stringstream stream;
  stream <<str;
  str.clear();
  stream >> v;
  vals.push_back(v);
  // Set up the array
  if (vals.empty()) {
    cout << "Empty array (?).";
    return;
  }
  size = vals.size();
  // Initialize
  initialize(size, sigma, phi);
  for (int i=0; i<size; i++) array[i] = vals.at(i);
  // Now test the solution
  solve(iters);
}

string Theory::print() {
  stringstream stream;
  stream << "{";
  for (int j=0; j<size; j++) {
    stream << "{" << j*scale << "," << array[j] << "}";
    if (j!=size-1) stream << ",";
  }
  stream << "}";
  string str;
  stream >> str;
  return str;
}

string Theory::printFreeLength() {
  stringstream stream;
  string str;
  stream << '{';
  for (int i=0; i<size; i++) {
    stream << "{" << i*scale << "," << lambda(i) << "}";
    if (i!=size-1) stream << ',';
  }
  stream << '}';
  stream >> str;
  return str;
}

string Theory::printFrequency() {
  stringstream stream;
  string str;
  stream << '{';
  for (int i=0; i<size; i++) {
    stream << '{';
    for (int j=0; j<size; j++) {
      stream << freq(i,j);
      if (j!=size-1) stream << ',';
    }
    stream << '}';
    if (i!=size-1) stream << ',';
  }
  stream << '}';
  stream >> str;
  return str;
}

string Theory::printFreqL() {
  stringstream stream;
  string str;
  stream << '{';
  for (int i=0; i<size; i++) {
    stream << "{" << i*scale << "," << freqL(i) << "}";
    if (i!=size-1) stream << ',';
  }
  stream << '}';
  stream >> str;
  return str;
}

string Theory::printFreqR() {
  stringstream stream;
  string str;
  stream << '{';
  for (int i=0; i<size; i++) {
    stream << "{" << i*scale << "," << freqR(i) << "}";
    if (i!=size-1) stream << ',';
  }
  stream << '}';
  stream >> str;
  return str;
}

string Theory::printPropagation() {
  double *deposit = new double[size];
  stringstream stream;
  string str;
  stream << '{';
  for (int i=0; i<size; i++) {
    stream << '{';
    getPropagation(i, freqL(i), freqR(i), deposit);
    //getPropagation(i, 0.5, 0.5, deposit);
    for (int j=0; j<size; j++) {
      stream << deposit[j];
      if (j!=size-1) stream << ',';
    }
    stream << '}';
    if (i!=size-1) stream << ',';
  }

  delete [] deposit;
  stream << '}';
  stream >> str;
  return str;
}

string Theory::printDArray() {
  stringstream stream;
  string str;
  stream << '{';
  for (int i=0; i<size; i++) {
    stream << "{" << i*scale << "," << d_array[i] << "}";
    if (i!=size-1) stream << ',';
  }
  stream << '}';
  stream >> str;
  return str;
}

void Theory::initialCondition() {
  // Uniform distribution

  //** NEED TO DO MORE HERE, CHANGE UNIFORM VALUE
  double uniform = phi/(PI*sqr(sigma));
  int start = sigma/scale;
  double max = 1/(2*sigma) / scale;
  array[start] = array[size-start-1] = max;
  for (int i=cutoff+1; i<size-cutoff-1; i++) array[i] = uniform;

  // ** FOR TESTING PURPOSES
  //for (int i=0; i<size; i++) array[i] = 0;
  //array[size-1] = 1/(2*sigma*scale); // Full

}

void Theory::getPropagation(int ri, double fL, double fR, double* deposit) {
  for (int i=0; i<size; i++) deposit[i] = 0; // Initialize [deposit] to 0's
  if (ri==0 || ri==size-1) return;
  // If you are hit from the left, you travel right and vice versa
  double PR = fL, PL = fR;
  // Rightwards propagation
  for (int rr=ri+1; rr<size-1; rr++) {
    double Pf = 1; // Probability of transmission
    // Integrate over the circumference of the disc
    for (int k=0; k<=propIters; k++) {
      double r2 = rr + height(k*sigma/propIters)/scale; // Plus
      Pf *= factor(r2);
    }
    Pf = pow(Pf, 1.0/(double)propIters);
    Pf *= (1-scale*drag); // Incorporate drag
    double diff = (1-Pf)*PR ; // This could still be wrong 
    if (diff<PR) {
      deposit[rr] += diff;
      PR -= diff;
    }
    else {
      deposit[rr] = PR;
      PR = 0;
      break;
    }
  }
  if (PR>0) deposit[size-1] += PR;
  
  // Leftwards propagation
  for (int rr=ri-1; 0<rr; rr--) {
    double Pf = 1;
    // Integrate over the circumference of the disc
    for (int k=0; k<=propIters; k++) {
      double r2 = rr - height(k*sigma/propIters)/scale; // Minus
      Pf *= factor(r2);
    }
    Pf = pow(Pf, 1.0/(double)propIters);
    Pf *= (1-scale*drag);
    double diff = (1-Pf)*PL;
    if (diff<PL) {
      deposit[rr] += diff;
      PL -= diff;
    }
    else {
      deposit[rr] = PL;
      PL = 0;
      break;
    }
  }
  if (PL>0) deposit[0] += PL;
}

double Theory::factor(int r) {
  if (r<0 || size<=r) return 0; // Particle sticks out of bounds
  double free = lambda(r);
  double factor = free>0 ? exp(-7*scale*(1/free-1)) : 0;
}

double Theory::vel(int i) {
  double R = 0.5*(size-1);
  return V0*(1-sqr((i-R)/R));
}

/// x is a BIN number, but the function can interpolate between bins
double Theory::gamma(double x) { 
  return fabs(scale*x*0.5) < sigma ? 2*sqrt(1-0.25*sqr(scale*x/sigma)) : 0;
}

double Theory::excL(double index) {
  int R = sigma/scale; // Radius, in bins
  double ex = 0;
  int start = max(0.,floor(index-R)), end = min((double)size-1.,ceil(index+R));
  // Integrate (we don't integrate with points finer then the bins so that we conserve particle number -> This was a problem before)
  for (int i=start; i<=end; i++) ex += array[i]*gamma(2*(index-i));
  ex *= (sigma*scale);
  return ex;
}

// Return -1 for inf
double Theory::freq(int r1, int r2) {
  if (fabs(r1-r2)>cutoff) return 0; // To far away
  double free = lambda(0.5*(r1+r2));
  double mult = 2*fabs(vel(r1)-vel(r2))*array[r1]*array[r2];
  if (mult==0) return 0;
  return free>0 ? mult/free : -1.;
}

// Frequency that an object at r is hit from the left
// Return -1 for inf
double Theory::freqL(int r) {
  double F=0;
  // Integrate
  for (int i=max(0,r-cutoff); i<r; i++) {
    double fr = freq(r,i);
    if (fr<0) return -1.;
    F += fr;
  }
  F*=scale; //
  return F;
}

// Frequency that an object at r is hit from the right
// Return -1 for inf
double Theory::freqR(int r) {
  double F=0;
  // Integrate
  for (int i=r+1; i<=min(size-1,r+cutoff); i++) {
    double fr = freq(r,i);
    if (fr<0) return -1.;
    F += fr;
  }
  F*=scale; //
  return F;
}

// Return -1 for inf
double Theory::freq(int r) {
  double F=0;
  for (int i=max(0,r-cutoff); i<min(size,r+cutoff); i++) {
    double fr = freq(r,i);
    if (fr<0) return -1;
    F += fr;
  }
  F*=scale; //
  return F;
}

double Theory::lambda(double r) {
  return 1-excL(r);
}

double Theory::at(double r) {
  int index = (int)r;
  double dx = r-index;
  return array[index]+(array[index+1]-array[index])*dx;
}

double Theory::height(double x) {
  return sqrt(sqr(sigma)-sqr(x));
}

double Theory::int_d_array() {
  double I = 0;
  for (int i=0; i<size; i++) I += fabs(d_array[i]);
  return I;
}
