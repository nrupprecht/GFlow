#ifndef __PRINTING_UTILITY_HPP__GFLOW__
#define __PRINTING_UTILITY_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  struct PrintingUtility {

    // Write an array that represents [elements]*[width] data to a .csv file
    static bool writeArrayDataToFile(RealType*, int, int, string);

    // [dirName] is the name of the directory that we should create our new directory of data in
    // [fileName] is the name of the new directory, and the files in that directory are called [fileName][#].csv
    static bool writeVectorToDirectory(vector<RealType*>&, const vector<int>&, int, string, const string);

    static bool writeVerletListToDirectory(const class VerletList&, const string);

    // Write a vector to a comma separated string
    static string toStrVec(RealType*);

    // Write an integer vector to a comma separated string
    static string toStrVec(int*);

  };

}
#endif // __PRINTING_UTILITY_HPP__GFLOW__