#ifndef THEORY_H
#define THEORY_H

#include "Utility.h"

class Theory {
 public:
  Theory(int, double, double);
  void initialize(int, double, double);

  // Solve function
  void solve(int);
  // Test function
  void test(string, double, int, double=1.); // Give an array as a string
  
  // Accessors
  double getRunTime() { return runTime; }

  string print();
  string printFreeLength(); 
  string printFrequency();
  string printFreqL();
  string printFreqR();
  string printPropagation();
  string printDArray();

  //private:
  // Helper functions
  void initialCondition();
  void getPropagation(int, double, double, double*);
  double factor(int r);
  double vel(int);
  double gamma(double); // Excluded length factor
  double excL(double);     // Excluded length (total)
  double freq(int,int); // Dependent frequency
  double freqR(int);    // Frequency that an object at r is hit from the right
  double freqL(int);    // Frequency that an object at r is hit from the left
  double freq(int);     // Total frequency
  double lambda(double);
  double at(double);
  double height(double);

  double int_d_array();

  double phi;
  double V0;
  double R0;   // Radius of the pipe
  double drag;
  double epsilon;
  double sigma; // Radius of the particles
  double scale;  // Length of a bin (in 'meters')
  int cutoff; // number of bins that corresponds to 2*radius
  int size;   // Number of bins
  double *array, *d_array;  

  int propIters;

  double runTime;
};

#endif
