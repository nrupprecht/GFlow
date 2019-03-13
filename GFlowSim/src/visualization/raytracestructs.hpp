#ifndef __RAY_TRACE_STRUCTS_HPP__GFLOW__
#define __RAY_TRACE_STRUCTS_HPP__GFLOW__

#include "../utility/utility.hpp"

namespace GFlowSimulation {

  /**
  *  \brief A data structure representing a ray, for ray tracing.
  *
  */
  struct Ray {
    //! \brief Constructor.
    Ray(const float*, const float*);

    //! \brief The origin point of the ray.
    float origin[3];

    //! \brief The orientation of the ray - this should be a normalized vector.
    float orientation[3];
  };

  /**
  *  \brief A data structure representing a sphere.
  *
  */
  struct Sphere {
    //! \brief Default constructor.
    Sphere();

    //! \brief Constructor.
    Sphere(const float*, const float);

    //! \brief Copy constructor.
    Sphere(const Sphere&);

    //! \brief Assignment operator.
    Sphere& operator=(const Sphere&);

    //! \brief Test whether a ray intersects a sphere. If so, return the intersection point.
    bool intersect(const Ray&, float*, float&, float&, const float) const;

    //! \brief A function that swaps two spheres.
    friend void swap(Sphere&, Sphere&);

    // --- Data

    //! \brief The center of the sphere.
    float center[3];

    //! \brief The radius of the sphere.
    float radius;
  };

}
#endif // __RAY_TRACE_STRUCTS_HPP__GFLOW__