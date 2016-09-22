#include "Theory.h"

Theory::Theory(int dim, double sigma, double phi) {
  epsilon = 0.005;
  radius = 1.;
  velocity = 1;
  Stat = dStat = Right = dRight = Left = dLeft = 0;
  initialize(dim, sigma, phi);
}

void Theory::initialize(int dim, double sigma, double phi) {
  
  this->dim = dim;
  this->sigma = sigma;
  this->phi = phi;

  scale = 2*radius/(double)dim;
  invScale = 1./scale;

  if (Stat) delete [] Stat;
  Stat = new double[dim];

  if (dStat) delete [] dStat;
  dStat = new double[dim];

  if (Right) delete [] Right;
  Right = new double[dim];

  if (dRight) delete [] dRight;
  dRight = new double[dim];

  if (Left) delete [] Left;
  Left = new double[dim];

  if (dLeft) delete [] dLeft;
  dLeft = new double[dim];

  // Initialize array
  double uniform = phi/(PI*sqr(sigma));
  int start = sigma*invScale;
  double max = invScale/(2*sigma);
  for (int i=0; i<dim; i++) Stat[i] = Right[i] = Left[i] = 0;
  Stat[start] = Stat[dim-start-1] = 0.9*max;
  int min = (int)(sqrt(3)*sigma*invScale);
  for (int i=start+min; i<dim-start-min; i++) Stat[i] = uniform;
}

void Theory::solve(int iters) {
  for (int it=0; it<iters; it++) update();
}

string Theory::print() {
  stringstream stream;
  stream << '{';
  for (int i=0; i<dim; i++) {
    stream << (Stat[i]+Right[i]+Left[i]);
    if (i!=dim-1) stream << ',';
  }
  stream << '}';
  string str;
  stream >> str;
  return str;
}

string Theory::printFreeLength() {
  stringstream stream;
  stream << '{';
  for (int i=0; i<dim; i++) {
    stream << lambda(i);
    if (i!=dim-1) stream << ',';
  }
  stream << '}';
  string str;
  stream >> str;
  return str;
}

string Theory::printFrequency() {
  stringstream stream;
  stream << '{';
  for (int i=0; i<dim; i++) {
    stream << freq(i);
    if (i!=dim-1) stream << ',';
  }
  stream << '}';
  string str;
  stream >> str;
  return str;
}

string Theory::printStat() {
  stringstream stream;
  stream << '{';
  for (int i=0; i<dim; i++) {
    stream << Stat[i];
    if (i!=dim-1) stream << ',';
  }
  stream << '}';
  string str;
  stream >> str;
  return str;
}

string Theory::printRight() {
  stringstream stream;
  stream << '{';
  for (int i=0; i<dim; i++) {
    stream << Right[i];
    if (i!=dim-1) stream << ',';
  }
  stream << '}';
  string str;
  stream >> str;
  return str;
}

string Theory::printLeft() {
  stringstream stream;
  stream << '{';
  for (int i=0; i<dim; i++) {
    stream << Left[i];
    if (i!=dim-1) stream << ',';
  }
  stream << '}';
  string str;
  stream >> str;
  return str;
}

string Theory::printDStat() {
  stringstream stream;
  stream << '{';
  for (int i=0; i<dim; i++) {
    stream << dStat[i];
    if (i!=dim-1) stream << ',';
  }
  stream << '}';
  string str;
  stream >> str;
  return str;
}

void Theory::update() {
  // Calculate updates
  for (int i=1; i<dim-1; i++) {
    dStat[i] = epsilon*(-freq(i)*Stat[i] + depR(i)*Right[i] + depL(i)*Left[i]);
    double DRDT = 0.5*invScale*(Right[i+1]+Right[i-1]-2*Right[i]);
    double dright = freqR(i)*Stat[i] - depR(i)*Right[i] - velocity*DRDT;
    dRight[i] = epsilon*dright;
    double DLDT = 0.5*invScale*(Left[i+1]+Left[i-1]-2*Left[i]);
    double dleft = freqL(i)*Stat[i] - depL(i)*Left[i] + velocity*DLDT;
    dLeft[i] = epsilon*dleft;
  }
  // Update
  for (int i=1;i<dim-1; i++) {
    Stat[i] += dStat[i];
    Right[i] += dRight[i];
    Left[i] += dLeft[i];
  }
}

double Theory::gamma(int bins) {
  return fabs(scale*bins*0.5) < sigma ? 2*sqrt(1-0.25*sqr(scale*bins/sigma)) : 0;
}

double Theory::excL(double index) {
  int R = sigma/scale; // Radius, in bins
  double ex = 0;
  int start = max(0.,floor(index-R)), end = min((double)dim-1.,ceil(index+R));
  // Integrate (we don't integrate with points finer then the bins so that we conserve particle number -> This was a problem before)
  for (int i=start; i<=end; i++) ex += (Stat[i]+Right[i]+Left[i])*gamma(2*(index-i));
  ex *= (sigma*scale);
  return ex;
}

double Theory::lambda(double r) {
  return 1-excL(r);
}

double Theory::vel(int bin) {
  double R = 0.5*(dim-1);
  return velocity*(1-sqr((bin-R)/R));
}

double Theory::freq(int r1, int r2) {
  if (abs(r1-r2)>2*sigma*invScale) return 0;
  return 2*fabs(vel(r1)-vel(r2))/lambda(0.5*(r1+r2));
}

double Theory::freqL(int r) {
  double F = 0;
  int end = min(dim-1, (int)(r+2*sigma*invScale));
  for (int i=r+1; i<end; i++) F += freq(i,r)*(Stat[i]+Left[i]+Right[i]);
  return scale*F;
}

double Theory::freqR(int r) {
  double F = 0;
  int start = max(0, (int)(r-2*sigma*invScale));
  for (int i=start; i<r; i++) F += freq(i,r)*(Stat[i]+Left[i]+Right[i]);
  return scale*F;
}

double Theory::freq(int r) {
  return freqL(r)+freqR(r);
}

double Theory::depR(int r) {
  return 0; //**
}

double Theory::depL(int r) {
  return 0; //**
}

