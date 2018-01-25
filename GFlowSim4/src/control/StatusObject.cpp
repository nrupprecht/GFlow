#include "StatusObject.hpp"

namespace GFlow {

  void StatusObject::writeMessage(std::string message) {
    messages.push_back(message);
  }

  void StatusObject::writeError(std::string error) {
    errors.push_back(error);
  }

  void StatusObject::writeArguments(std::string args) {
    arguments = args;
  }

}
