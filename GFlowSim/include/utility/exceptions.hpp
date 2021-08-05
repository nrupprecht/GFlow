#ifndef __EXCEPTIONS_HPP__GFLOW__
#define __EXCEPTIONS_HPP__GFLOW__

namespace GFlowSimulation {

/**
*  @brief The base class for GFlow specific exceptions.
*
*  The exception contains a string that we use to write a
*  message containing more information about the error that
*  occured.
*/
struct Exception : public std::exception {
    //! Default constructor.
    Exception() {};

    //! Message constructor.
    Exception(const std::string &m) : message(m) {};

    //! \brief Override the what function.
    const char *what() const noexcept override {
        return message.c_str();
    }

    //! \brief The error message.
    std::string message;
};

/**
*  \brief Exception for encountering an unexpected null pointer.
*
*  This exception is thrown whenever we encounter a null pointer
*  when we really shouldn't. E.g. if the object should have always
*  been initialized before getting to this point in the simulation.
*/
struct UnexpectedNullPointer : public Exception {
    //! Default constructor.
    UnexpectedNullPointer() {};

    //! Message constructor.
    UnexpectedNullPointer(const std::string &m) : Exception(m) {};
};

/**
*  @brief Exception for encountering bad or mismatched dimensions.
*
*  This exception is thrown whenever we encounter an unexpected
*  difference between the dimensionality of objects.\n
*  There is a separate exception for encountering or creating a
*  zero dimensional object, which should never occur.
*/
struct BadDimension : public Exception {
    //! Default constructor.
    BadDimension() {};

    //! Message constructor.
    BadDimension(const std::string &m) : Exception(m) {};
};

/**
*  @brief Exception for encountering a zero dimensional object.
*
*  This exception is thrown whenever we encounter a zero dimensional
*  object. There are multiple template<int> objects that should never
*  be zero dimensional. This error is thrown if one such object is
*  zero dimensional.
*/
struct ZeroDimension : public Exception {
    //! Default constructor.
    ZeroDimension() {};

    //! Message constructor.
    ZeroDimension(const std::string &m) : Exception(m) {};
};

/**
*  @brief Exception for encountering particles of types that shouldn't
*  exist.
*
*  This exception is thrown whenever we encounter a particle of type
*  greater than ntypes - 1, since only particles of type -1 (not a
*  particle), 0, 1, ... , (ntypes-1) should exist.
*/
struct ParticleTypeError : public Exception {
    //! Default constructor.
    ParticleTypeError() {};

    //! Message constructor.
    ParticleTypeError(const std::string &m) : Exception(m) {};
};

/**
*  @brief Exception for calling an unimplemented function.
*
*  This exception is thrown whenever we try to call a function
*  that should not be called because it is unimplemented, and
*  needs to be implemented for proper functionality.
*/
struct Unimplemented : public Exception {
    //! \brief Default constructor.
    Unimplemented() {};

    //! \brief Message constructor.
    Unimplemented(const std::string &m) : Exception(m) {};
};

/**
*  @brief Exception for when a dataF or dataI array is missing.
*/
struct DataEntryNotFound : public Exception {
    //! Default constructor.
    DataEntryNotFound() {};

    //! Message constructor.
    DataEntryNotFound(const std::string &m) : Exception(m) {};
};

/**
*  @brief Exception for when a value is a Nan.
*/
struct NanValue : public Exception {
    //! Default constructor.
    NanValue() {};

    //! Message constructor.
    NanValue(const std::string &m) : Exception(m) {};
};

/**
*  @brief Exception for when bounds are bad.
*/
struct BadBounds : public Exception {
    //! Default constructor.
    BadBounds() {};

    //! Message constructor.
    BadBounds(const std::string &m) : Exception(m) {};
};

}
#endif // __EXCEPTIONS_HPP__GFLOW__