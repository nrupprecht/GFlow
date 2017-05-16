#include "FileParser.hpp"

namespace GFlow {
  
  bool FileParser::parser(string filename) {
    ifstream fin(filename);
    if (fin.fail()) {
      return false;
    }

    char c;
    fin.get(c);
    while (!fin.eof()) {

      
      // (comment) // 
      if (c=='/') fin.getline();
      
      // wrapX = [true/false]
      
      // wrapY = [true/false]
      
      // wall - [lx] [ly] [rx] [ry]
      
      // region - 
      // [left] [right] [bottom] [top]
      // (options) N=[number] | disp=[dispersion] | disp_type=[dispersion type] | it=[interaction] | ...
      // end
      
      fin.get(c);
    }

    return true;
  }
  
}
