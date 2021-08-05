#include "utility/argparse2.hpp"

ArgParse::ArgParse(int ac, char **av)
    : argc(ac), argv(av) {};

// Parse the command line arguments for information
void ArgParse::parse() {
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') { // Our cue that this is a token
      string tok, val;
      int j = 1;
      char c = argv[i][j];
      while (c != '\0' && c != '=') {
        tok.push_back(c); // Record the token
        j++;
        c = argv[i][j];
      }
      if (c != '\0') {
        j++;
        c = argv[i][j];
      }
      while (c != '\0') { // Get the value
        val.push_back(c);
        j++;
        c = argv[i][j];
      }
      tlist.push_back(pair<std::string, std::string>(tok, val));
      checked.push_back(false);
    }
    else {
      throw IllegalToken(argv[i][0]);
    }
  }
}

void ArgParse::specifyHelpCommand(const std::string& help_command) {

}

std::string ArgParse::printHelpMessage() {
  // TODO: Stub.
  return "";
}

std::string ArgParse::checkForMultiples() {
  // TODO: Stub.
  return "";
}

bool ArgParse::check() {
  // @todo Stub.
  return true;
}