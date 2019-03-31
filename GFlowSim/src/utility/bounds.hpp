#ifndef __BOUNDS_HPP__GFLOW__
#define __BOUNDS_HPP__GFLOW__

#include "utility.hpp"

namespace GFlowSimulation {

  /*
  *  @struct Bounds
  *
  *  Standard object for storing rectangular bounds of arbitrary dimensionality.
  *  NOTE: For safety, we do not include a move equals operator or move constructor.
  *
  */
  struct Bounds {
    //! \brief Dimension setting constructor.
    Bounds(const int);

    //! \brief Full setting constructor.
    Bounds(const RealType*, const RealType*, const int);

    //! \brief Copy constructor
    Bounds(const Bounds&);

    //! \brief Destructor.
    ~Bounds();

    //! \brief Copy equals.
    Bounds& operator=(const Bounds&);

    //! \brief Boolean equals operator
    bool operator==(const Bounds&) const;

    //! \brief Boolean not equals operator
    bool operator!=(const Bounds&) const;

    //! \brief Check whether the bounds contains a point.
    bool contains(RealType*) const;

    //! \brief Get the widths in various dimensions of the bounds.
    RealType wd(int) const;

    //! \brief Get the volume of the bounds.
    RealType vol() const;

    //! \brief Get the dimensionality of the bounds.
    int dims() const;

    //! \brief Set the input vector to be the center of the bounds.
    void center(RealType*) const;

    //! \brief Set the input vector to be a uniformly random position from within the bounds. 
    void randomPoint(RealType*) const;

    //! \brief Find the maximum bounds width.
    friend RealType max_width(const Bounds&);

    //! \brief Find the minimum bounds width.
    friend RealType min_width(const Bounds&);

    // --- Data

    //! \brief The min and max.
    RealType *min, *max;

    //! \brief The dimensionality of this bounds object
    int dimensions;
  };

}

#endif // __BOUNDS_HPP__GFLOW__
