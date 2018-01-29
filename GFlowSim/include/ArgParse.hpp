#ifndef ARG_PARSE_H
#define ARG_PARSE_H

#include <string>
using std::string;

#include <sstream>
using std::stringstream;

#include <vector>
using std::pair;
using std::vector;

#include <iostream>
using std::cout;
using std::endl;

/// Command line argument parsing class
class ArgParse {
 public:
  // Default constructor
  ArgParse() : argc(0), argv(0) {};
  
  // Constructor that sets argc, argv
  ArgParse(int argc, char** argv) {
    this->argc = argc;
    this->argv = argv;
    parse();
    throw true; 
  }
  
  // Set argc, argv
  void set(int ac, char** av) { 
    argc = ac; 
    argv = av; 
    // Parse the arguments
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
	checked.push_back(false);
      }
      else throw IllegalToken(argv[i][0]);
    }
  }

  // Returns a pair of empty strings if not found
  pair<string,string> find(string token) {
    pair<string,string> opt("","");
    int i=0;
    for (auto tpair : tlist) {
      if (tpair.first==token) {
	opt = tpair;
	// We have looked for this token
	checked.at(i) = true;
      }
      ++i;
    }
    return opt;
  }

  // The token to look for and a variable in which to store the value it was given
  template<typename T> void get(const string& token, T& var) {
    // Find the {token, value} pair
    pair<string,string> opt = find(token);
    // If the token was found
    if (!opt.first.empty()) {
      stringstream stream;
      // If no argument was provided, use 1 (true) as the default argument
      if (opt.second.empty()) stream << "1";
      // Otherwise, read in the argument. Stringstream converts it to the correct type
      else stream << opt.second;
      stream >> var;
    }
  }

  // Checks if any command line arguments were not checked by the parser. If so, then they are assumed to be illegal arguments, and an error will be thrown
  void check() {
    for (int i=0; i<checked.size(); ++i) {
      if (!checked.at(i)) throw UncheckedToken("Tok=["+tlist.at(i).first + "]: Val=[" + tlist.at(i).second+"]");
    }
  }

  /// Exception classes
  class IllegalToken {
  public:
    IllegalToken(char d) : c(d) {};
    char c;
  };

  /// Exception class
  class UncheckedToken {
  public:
    UncheckedToken(string s) : token(s) {};
    string token;
  };

private:
  // Number of arguments and string array
  int argc;
  char** argv;
  
  // List of found tokens and their corresponding value
  vector<pair<string,string> > tlist;

  // Tokens we have looked for in the argv
  vector<bool> checked;

};

#endif
