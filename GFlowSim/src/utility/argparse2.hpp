#ifndef __ARGPARSE_2_HPP__
#define __ARGPARSE_2_HPP__

#include <map>
using std::multimap;

#include <string>
using std::string;

#include <vector>
using std::vector;
using std::pair;

class ArgParse {
public:
  //! @brief Default constructor.
  ArgParse() : argc(0), argv(nullptr) {};

  //! @brief Setting constructor.
  ArgParse(int, char**);

  //! @brief Parse the command line arguments.
  void parse();

  //! @brief Specify which token activates the help message.
  void specifyHelpCommand(string);

  //! @brief Print out a help message.
  string printHelpMessage();

  //! @brief Add a token, and what variable to fill with the value.
  //!
  //! Also can add an explanation string
  template<typename T> void get(const string token, T& var, string explanation="") {
    // Record the explanation, if there is one.
    if (!explanation.empty()) explanations.push_back(pair<string, string>(token, explanation));

    auto it = cmdlineargs.find(token);
  } 

  //! @brief Check if any arguments have been given multiple values. 
  //!
  //! Print out a message string if multiples are found.
  string checkForMultiples();

  //! @brief Check if any arguments were in the command line that we didn't check for.
  bool check();

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

  //! @brief A list of tokens and values.
  multimap<string, string> cmdlineargs;

  //! @brief Tokens, and their explanation (for the help message).
  vector<pair<string, string> > explanations;

  //! @brief Which token activates the help message.
  string helpToken = "";

  //! @brief The number of arguments.
  int argc;
  //! @brief The arguments.
  char **argv;
}


#endif // __ARGPARSE_2_HPP__