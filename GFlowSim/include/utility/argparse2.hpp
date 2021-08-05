#ifndef __ARGPARSE_2_HPP__
#define __ARGPARSE_2_HPP__

#include <map>
#include <string>
#include <vector>

#include "utility/exceptions.hpp"

class ArgParse {
 public:
  //! \brief Default constructor.
  ArgParse() = default;

  //! @brief Setting constructor.
  ArgParse(int, char **);

  //! \brief Parse the command line arguments.
  void parse();

  //! \brief Specify which token activates the help message.
  void specifyHelpCommand(const std::string &command);

  //! \brief Print out a help message.
  std::string printHelpMessage();

  //! \brief Add a token, and what variable to fill with the value.
  //!
  //! Also can add an explanation string
  template<typename T>
  void get(const std::string &token, T &var, const std::string &explanation = "") {
    // Record the explanation, if there is one.
    if (!explanation.empty()) {
      explanations.emplace_back(token, explanation);
    }

    auto it = cmdlineargs.find(token);
  }

  //! @brief Check if any arguments have been given multiple values. 
  //!
  //! Print out a message string if multiples are found.
  std::string checkForMultiples();

  //! @brief Check if any arguments were in the command line that we didn't check for.
  bool check();

  /// Exception classes
  class IllegalToken : public std::exception {
   public:
    explicit IllegalToken(char d)
        : c(d) {};
    char c;

    const char *what() const noexcept override {
      return (std::string("Illegal token: ") + c).c_str();
    }
  };

  /// Exception class
  class UncheckedToken : public GFlowSimulation::Exception {
   public:
    explicit UncheckedToken(const std::string &s)
        : Exception(s) {};
  };

 private:

  //! \brief A list of tokens and values.
  std::multimap<std::string, std::string> cmdlineargs{};

  //! \brief Tokens, and their explanation (for the help message).
  std::vector<std::pair<std::string, std::string>> explanations{};

  //! \brief Which token activates the help message.
  std::string helpToken{};

  //! \brief The number of arguments.
  int argc = 0;
  //! \brief The arguments.
  char **argv = nullptr;
};

#endif // __ARGPARSE_2_HPP__