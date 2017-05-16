#ifndef __DATA_CREATOR_HPP__
#define __DATA_CREATOR_HPP__

// Includes
#include "../../include/Utility.hpp"

namespace GFlow {

  // Forward declaration to region
  struct Region;

  // Data creation class
  template<typename T> struct DataCreator {
    virtual vector<T> makeValues(Region*);
  };

  class Fixed_Phi_Uniform_Radii : public DataCreator<RealValue> {
  public:
    Fixed_Phi_Uniform_Radii(RealType ph, RealType sig, RealType disp) : phi(ph), sigma(sig), dispersion(disp) {};

    virtual vector<RealType> makeValues(Region*);

  private:
    RealType phi, sigma, dispersion;
  };

}
#endif // __DATA_CREATOR_HPP__
