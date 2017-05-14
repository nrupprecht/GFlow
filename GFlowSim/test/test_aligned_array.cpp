/*
 * @author Nathaniel Rupprecht
 * @created May 13, 2017
 * Testing the aligned_array class
 *
 */

// Includes
#include "../include/Utility.hpp"
#include "../include/aligned_array.hpp"
using GFlow::aligned_array;

int main(int argc, char** argv) {
  // Create an array
  aligned_array<double> array;
  // Reserve size of 100
  array.reserve(100);
  // Print initial values
  for (int i=0; i<100; ++i) cout << array.at(i) << " ";
  cout << endl;
  // Assign new values
  for (int i=0;i<100; ++i) array[i] = i*i;
  cout << endl;
  // Print final values
  for (int i=0;i<100; ++i) cout << array.at(i) << " ";

  return 0;
}
