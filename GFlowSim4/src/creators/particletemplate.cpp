#include "particletemplate.hpp"

namespace GFlowSimulation {

  ParticleTemplate::ParticleTemplate() : type_engine(nullptr), radius_engine(nullptr), mass_engine(nullptr), velocity_engine(nullptr),
    type_string(""), radius_string(""), mass_string(""), velocity_string("") {};

  //! @brief Destructor.
  ParticleTemplate::~ParticleTemplate() {
    if (type_engine)     delete type_engine;
    if (radius_engine)   delete radius_engine;
    if (mass_engine)     delete mass_engine;
    if (velocity_engine) delete velocity_engine;
  }

  ParticleTemplate::ParticleTemplate(const ParticleTemplate& p) {
    type_engine     = p.type_engine ? p.type_engine->copy() : nullptr;
    radius_engine   = p.radius_engine ? p.radius_engine->copy() : nullptr;
    mass_engine     = p.mass_engine ? p.mass_engine->copy() : nullptr;
    velocity_engine = p.velocity_engine ? p.velocity_engine->copy() : nullptr;
    type_string = p.type_string;
    radius_string = p.radius_string;
    mass_string = p.mass_string;
    velocity_string = p.velocity_string;
  }

  ParticleTemplate::ParticleTemplate(ParticleTemplate&& p) {
    type_engine     = p.type_engine;
    radius_engine   = p.radius_engine;
    mass_engine     = p.mass_engine;
    velocity_engine = p.velocity_engine;
    p.type_engine = nullptr;
    p.radius_engine = nullptr;
    p.mass_engine = nullptr;
    p.velocity_engine = nullptr;
    type_string = p.type_string;
    radius_string = p.radius_string;
    mass_string = p.mass_string;
    velocity_string = p.velocity_string;
  }

  ParticleTemplate& ParticleTemplate::operator=(const ParticleTemplate& p) {
    type_engine     = p.type_engine->copy();
    radius_engine   = p.radius_engine->copy();
    mass_engine     = p.mass_engine->copy();
    velocity_engine = p.velocity_engine->copy();
    type_string = p.type_string;
    radius_string = p.radius_string;
    mass_string = p.mass_string;
    velocity_string = p.velocity_string;
    return *this;
  }

  void ParticleTemplate::createParticle(RealType *X, RealType& radius, RealType &im, int& type, int n) {
    if (type_engine==nullptr || radius_engine==nullptr || mass_engine==nullptr) {
      cout << type_engine << " " << radius_engine << " " << mass_engine << endl;
      throw NullEngine();
    }
    type = static_cast<int>(type_engine->generate());
    radius = radius_engine->generate();
    double m = mass_engine->generate();

    cout << m << endl;

    if (mass_string=="Density") m = sphere_volume(radius)*m;

    cout << mass_string << " " << m << " " << sphere_volume(radius) << " " << radius << endl;

    im = m>0 ? 1./m : 0;
  }

}