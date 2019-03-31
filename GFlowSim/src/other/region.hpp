#ifndef __REGION_HPP__GFLOW__
#define __REGION_HPP__GFLOW__

#include "../utility/bounds.hpp"

namespace GFlowSimulation {

  /** 
  *  \brief The base class for spatial region objects.
  *
  */
  class Region {
  public:
    //! \brief Constructor. Sets the dimensionality.
    Region(int d) : dimensions(d) {};

    //! \brief Doesn't do anything, but keeps warnings from occuring.
    virtual ~Region() {};

    //! \brief Pick a position at uniform random from within the bounds.
    virtual void pick_position(RealType*) = 0;

    //! \brief Return whether a point falls within the bounds.
    virtual bool contains(RealType*) = 0;

    //! \brief Return the volume of the bounds.
    virtual RealType vol() = 0;

    //! \brief Returns bounds that serve as a bounding box for the region.
    virtual Bounds get_bounding_box() = 0;

    //! \brief Return the dimensionality of the region.
    int get_dimensions() { return dimensions; }

  protected:
    //! \brief The dimensionality of the bounds.
    const int dimensions;
  };


  class RectangularRegion : public Region {
  public:
    //! \brief Dimension setting constructor.
    RectangularRegion(int);

    //! \brief Dimension and bounds setting constructor.
    RectangularRegion(int, const Bounds&);

    //! \brief Pick a position at uniform random from within the bounds.
    virtual void pick_position(RealType*) override;

    //! \brief Return whether a point falls within the bounds.
    virtual bool contains(RealType*) override;

    //! \brief Return the volume of the bounds.
    virtual RealType vol() override;

    //! \brief Returns bounds that serve as a bounding box for the region.
    virtual Bounds get_bounding_box() override;

    // --- Accessors

    //! \brief Get (by reference) the bounds min values.
    RealType& min(int);

    //! \brief Get (by reference) the bounds max values.
    RealType& max(int);
    
    // --- Mutators

    //! \brief Sets the bounds of the rectangular region.
    void setBounds(const Bounds&);

  private:
    //! \brief A rectangular region can be represented by a Bounds object.
    Bounds bounds;
  };

  class SphericalRegion : public Region {
  public:
    //! \brief Dimension setting constructor.
    SphericalRegion(int);

    //! \brief Dimension, center, and radius setting constructor.
    SphericalRegion(int, const RealType*, const RealType);

    //! \brief Destructor.
    virtual ~SphericalRegion();

    //! \brief Pick a position at uniform random from within the bounds.
    virtual void pick_position(RealType*) override;

    //! \brief Return whether a point falls within the bounds.
    virtual bool contains(RealType*) override;

    //! \brief Return the volume of the bounds.
    virtual RealType vol() override;

    //! \brief Returns bounds that serve as a bounding box for the region.
    virtual Bounds get_bounding_box() override;

    // --- Accessors

    //! \brief Access a component of the region's center vector.
    RealType &center(int);

    //! \brief Access the radius of the region.
    RealType &radius();

    // --- Mutators

    //! \brief Set the center of the region.
    void setCenter(const RealType*);

    //! \brief Set the radius of the region.
    void setRadius(const RealType);

  private:
    //! \brief The center of the region.
    RealType *_center;

    //! \brief The radius of the region.
    RealType _radius;
  };

}
#endif // __REGION_HPP__GFLOW__