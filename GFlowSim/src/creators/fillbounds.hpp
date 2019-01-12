#ifndef __FILL_BOUNDS_HPP__GFLOW__
#define __FILL_BOUNDS_HPP__GFLOW__

#include "../base/creator.hpp"

namespace GFlowSimulation {

  //! \brief The base class for fill bounds.
  //!
  //! Fill bounds are used to specify the bounds of a fill-area for the file parse creator.
  //! \see FileParseCreator
  struct FillBounds {
    FillBounds(int dim) : bnd_dimensions(dim) {};

    //! \brief Virtual destructor. 
    //!
    //! Doesn't do anything, but keeps warnings from happening.
    virtual ~FillBounds() {};

    //! \brief Returns the volume of the bounds
    virtual double vol()=0;
    //! \brief Return a random position within the bounds
    virtual void pick_position(RealType*)=0;
    //! \brief Get rectangular bounds that enclose the bounds
    virtual Bounds getBounds()=0;

    //! \brief Number of dimensions the bounds takes up.
    //!
    //! We include this so we can make lower dimensional fill areas - walls, lines, etc
    const int bnd_dimensions;
  };

  //! \brief Rectangular fill bounds.
  struct RectangularBounds : public FillBounds {
    RectangularBounds(int);

    RectangularBounds(Bounds, int);

    virtual ~RectangularBounds();

    virtual double vol() override;

    virtual void pick_position(RealType *x) override;

    Bounds getBounds() override;

    float *min, *max;
  };

  //! \brief Spherical fill bounds.
  struct SphericalBounds : public FillBounds {
    SphericalBounds(int);

    virtual ~SphericalBounds();

    virtual double vol() override;

    virtual void pick_position(RealType *x) override;

    Bounds getBounds() override;

    float *center;
    float radius;
  };

}
#endif // __FILL_BOUNDS_HPP__GFLOW__
