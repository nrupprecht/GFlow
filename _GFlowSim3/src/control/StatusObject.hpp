#ifndef __STATUS_OBJECT_HPP__
#define __STATUS_OBJECT_HPP__

// Includes
#include <string>
#include <vector>

namespace GFlow {
  
  /*
   * @class StatusObject
   *
   */
  class StatusObject{
  public:
    // Accessors
    std::vector<std::string> getMessages() const { return messages; }
    std::vector<std::string> getErrors()   const { return errors; }
    std::string getArguments()        const { return arguments; }

    // Mutators
    void writeMessage(std::string);
    void writeError(std::string);
    void writeArguments(std::string);

  private:
    std::vector<std::string> messages;
    std::vector<std::string> errors;
    std::string arguments;
  };

};
#endif // __STATUS_OBJECT_HPP__
