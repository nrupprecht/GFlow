#ifndef THEORY_H
#define THEORY_H

#include "Utility.h"

class Theory {
 public:
  Theory(int, double, double);
  void initialize(int, double, double);

  void solve(int);     // Run for some number of iterations

  double getRunTime() { return 0; }

  // Printing functions
  string print();
  string printFreeLength();
  string printFrequency();
  string printStat();
  string printRight();
  string printLeft();
  string printDStat();

 private:
  
  void update();     // Update the distributions for a time step
  
  // Helper functions
  double gamma(int); // Geometry factor with the argument in NUMBER OF BINS (not length in 'meters')
  double excL(double);   // Excluded length
  double lambda(double); // Free length ( 1 - excluded length )
  double vel(int);   // Velocity as a function of BINS (not position in 'meters')
  double excL(int);  // Excluded length
  double freq(int, int); // Collision factor between bins
  double freqL(int);     // Leftwards collision factor
  double freqR(int);     // Rightwards collision factor
  double freq(int);      // Total collision factor (left + right collision factors)
  double depR(int);      // Rightwards deposition probability
  double depL(int);      // Leftwards deposition probability
  
  // Dimensions
  int dim;             // Number of bins to use
  double sigma;        //  Disc radius
  double phi;          //
  double epsilon;      // Iteration time step
  double radius;       // Radius of the pipe
  double scale;        // The length of a single bin
  double invScale;     // The inverse of binLength
  double velocity;     // The velocity of moving balls

  // Distributions and change in distributions
  double *Stat, *Right, *Left;
  double *dStat, *dRight, *dLeft;
};

#endif // THEORY_H
