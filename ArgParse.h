#ifndef ARG_PARSE_H
#define ARG_PARSE_H

#include <string>
using std::string;

#include <sstream>
using std::stringstream;

#include <vector>
using std::pair;
using std::vector;

/// Command line argument parsing class
class ArgParse {
public:
  ArgParse(int argc, char** argv) {
    this->argc = argc;
    this->argv = argv;
    parse();
  }
  
  // Parse the command line arguments for information
  void parse() {
    for (int i=1; i<argc; i++) {
      if (argv[i][0]=='-') { // Our cue that this is a token
	string tok, val;
	int j=1;
	char c = argv[i][j];
	while(c!='\0' && c!='=') {
	  tok.push_back(c); // Record the token
	  j++;
	  c = argv[i][j];
	}
	if (c!='\0') {
	  j++;
	  c = argv[i][j];
	}
	while(c!='\0') { // Get the value
	  val.push_back(c);
	  j++;
	  c = argv[i][j];
	}
	tlist.push_back(pair<string,string>(tok,val));
      }
      else throw IllegalToken();
    }
  }

  // Returns a pair of empty strings if not found
  pair<string,string> find(string token) {
    pair<string,string> opt("","");
    for (auto tpair : tlist) {
      if (tpair.first==token) {
	opt = tpair;
      }
    }
    return opt;
  }

  template<typename T> void get(const string& token, T& var) {
    pair<string,string> opt = find(token);
    if (!opt.first.empty()) {
      stringstream stream;
      if (opt.second.empty()) stream << "1"; // Default argument
      else stream << opt.second;
      stream >> var;
    }
  }

  /// Exception classes
  class IllegalToken {};

private:
  int argc;
  char** argv;
  
  vector<pair<string,string> > tlist; // List of found tokens and their corresponding value

};

#endif
