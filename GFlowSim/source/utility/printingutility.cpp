#include "utility/printingutility.hpp"

namespace GFlowSimulation {

bool writeArrayDataToFile(RealType *array,
                          int elements,
                          int width,
                          const std::string &fileName) {
  // Open a file stream
  ofstream fout(fileName);
  if (fout.fail()) {
    return false;
  }

  // Write the data
  for (int n = 0; n < elements; ++n) {
    // Print out:  x[0], x[1], ...  \n
    for (int d = 0; d < width; ++d) {
      fout << array[width * n + d];
      if (d != width - 1) {
        fout << ",";
      }
    }
    fout << "\n";
  }

  // Return success
  return true;
}

//! \param dirName is the name of the directory that we should create our new directory of data in
//! \param fileName is the name of the new directory, and the files in that directory are called [fileName][#].csv
bool writeVectorToDirectory(std::vector<RealType *> &record,
                            const std::vector<int> &elements,
                            int width,
                            std::string dirName,
                            const std::string &fileName) {
  // Create the directory
  if (*dirName.rbegin() == '/') { // Make sure there is a /
    dirName += fileName + "/";
  }
  else {
    dirName += ("/" + fileName + "/");
  }
  // --> Dir name is now the name of the directory we are creating
  // Create a directory for all the data
  mkdir(dirName.c_str(), 0777);

  // Create the data files in the directory
  for (uint iter = 0; iter < record.size(); ++iter) {
    // The name of the file for this data
    string subfile = dirName + fileName + toStr(iter) + ".csv";
    // Create the file [fileName] (inside directory [dirName]) containing the data from record.at(iter)
    if (!writeArrayDataToFile(record.at(iter), elements.at(iter), width, subfile)) {
      return false;
    }
  }

  // Return success
  return true;
}

std::string toStrVec(const float *x, int length) {
  std::string str;
  for (int d = 0; d < length; ++d) {
    str += toStr(x[d]);
    if (d != length - 1) {
      str += ',';
    }
  }
  return str;
}

std::string toStrVec(const double *x, int length) {
  std::string str;
  for (int d = 0; d < length; ++d) {
    str += toStr(x[d]);
    if (d != length - 1) {
      str += ',';
    }
  }
  return str;
}

std::string toStrVec(const int *x, int length) {
  std::string str;
  for (int d = 0; d < length; ++d) {
    str += toStr(x[d]);
    if (d != length - 1) {
      str += ',';
    }
  }
  return str;
}

}