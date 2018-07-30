#ifndef __EXCEPTIONS_HPP__GFLOW__
#define __EXCEPTIONS_HPP__GFLOW__

namespace GFlowSimulation {

  struct Exception {
    // Default constructor
    Exception() {};
    
    // Message constructor
    Exception(const string& m) : message(m) {};

    string message;
  };

  struct UnexpectedNullPointer : public Exception {
    UnexpectedNullPointer() {};
    UnexpectedNullPointer(const string& m) : Exception(m) {};
  };

  struct BadDimension : public Exception {
    BadDimension() {};
    BadDimension(const string& m) : Exception(m) {};
  };

  struct ParticleTypeError : public Exception {
    ParticleTypeError() {};
    ParticleTypeError(const string& m) : Exception(m) {};
  };

  struct Unimplemented : public Exception {
    Unimplemented() {};
    Unimplemented(const string& m) : Exception(m) {};
  };

}
#endif // __EXCEPTIONS_HPP__GFLOW__