#include "printingutility.hpp"
// Other files
#include "verletlist.hpp"

namespace GFlowSimulation {

  bool PrintingUtility::writeArrayDataToFile(RealType *array, int elements, int width, string fileName) {
    // Open a file stream
    ofstream fout(fileName);
    if (fout.fail()) return false;

    // Write the data
    for (int n=0; n<elements; ++n) {
      // Print out:  x[0], x[1], ...  \n
      for (int d=0; d<width; ++d) {
        fout << array[width*n+d];
        if (d!=width-1) fout << ",";
      }
      fout << "\n";
    }

    // Return success
    return true;
  }

  // [dirName] is the name of the directory that we should create our new directory of data in
  // [fileName] is the name of the new directory, and the files in that directory are called [fileName][#].csv
  bool PrintingUtility::writeVectorToDirectory(vector<RealType*>& record, const vector<int>& elements, 
    int width, string dirName, const string fileName) 
  {
    // Create the directory
    if (*dirName.rbegin()=='/') // Make sure there is a /
      dirName += fileName+"/";
    else 
      dirName += ("/"+fileName+"/");
    // --> Dir name is now the name of the directory we are creating

    // Create a directory for all the data
    mkdir(dirName.c_str(), 0777);

    // Create the data files in the directory
    for (uint iter=0; iter<record.size(); ++iter) {
      // The name of the file for this data
      string subfile = dirName+fileName+toStr(iter)+".csv";
      // Create the file [fileName] (inside directory [dirName]) containing the data from record.at(iter)
      if (!writeArrayDataToFile(record.at(iter), elements.at(iter), width, subfile)) return false;
    }

    // Return success
    return true;
  }

  bool PrintingUtility::writeVerletListToDirectory(const VerletList &vl, const string fileName) {
    // Open file stream
    ofstream fout(fileName);
    if (fout.fail()) return false;

    // Get the data
    const int hsize = vl.vlHSize(), vsize = vl.vlSize();
    const int *heads = vl.getHeads(), *verlet = vl.getVerlet();
    // Write the data
    for (int h=0; h<hsize-1; ++h) {
      fout << verlet[heads[h]] << ",";
      for (int p=heads[h]+1; p<heads[h+1]; ++p) {
        fout << verlet[p];
        if (p!=heads[h+1]-1) fout << ",";
      }
      fout << endl;
    }
    // Last one
    fout << verlet[heads[hsize-1]] << ",";
    for (int p=heads[hsize-1]+1; p<vsize; ++p) {
      fout << verlet[p];
      if (p!=vsize-1) fout << ",";
    }

    // Return success
    return true;
  }

  string PrintingUtility::toStrVec(RealType *x) {
    string str;
    for (int d=0; d<DIMENSIONS; ++d) {
      str += toStr(x[d]);
      if (d!=DIMENSIONS-1) str += ',';
    }
    return str;
  }

  string PrintingUtility::toStrVec(int *x) {
    string str;
    for (int d=0; d<DIMENSIONS; ++d) {
      str += toStr(x[d]);
      if (d!=DIMENSIONS-1) str += ',';
    }
    return str;
  }

}
