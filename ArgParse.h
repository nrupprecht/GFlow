#ifndef ARG_PARSE_H
#define ARG_PARSE_H

#include <string>
using std::string;

#include <sstream>
using std::stringstream;

/// Command line argument parsing class
class ArgParse{
public:
  ArgParse(int argc, char** argv) {
    this->argc = argc;
    this->argv = argv;
    parse();
  }
  
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

private:
  int argc;
  char** argv;
  
  vector<pair<string,string> > tlist; // List of found tokens and their corresponding value

};

#endif
