#ifndef __SIM_DATA_HPP__
#define __SIM_DATA_HPP__

#include "SimDataBase.hpp"

namespace GFlow {

  class SimData : public SimDataBase<SimData> {
  public:
    SimData(const Bounds&, const Bounds&);

    // Make the base class a friend
    friend class SimDataBase<SimData>;
    
    // Pointer access to arrays
    RealType* getPxPtr() { return px.getPtr(); }
    RealType* getPyPtr() { return py.getPtr(); }
    RealType* getVxPtr() { return vx.getPtr(); }
    RealType* getVyPtr() { return vy.getPtr(); }
    RealType* getFxPtr() { return fx.getPtr(); }
    RealType* getFyPtr() { return fy.getPtr(); }
    RealType* getThPtr() { return th.getPtr(); }
    RealType* getOmPtr() { return om.getPtr(); }
    RealType* getSgPtr() { return sg.getPtr(); }
    RealType* getTqPtr() { return tq.getPtr(); }
    RealType* getImPtr() { return im.getPtr(); }
    RealType* getIiPtr() { return iI.getPtr(); }
    RealType* getRpPtr() { return rp.getPtr(); }
    RealType* getDsPtr() { return ds.getPtr(); }
    RealType* getCfPtr() { return cf.getPtr(); }
    int*      getItPtr() { return it.getPtr(); }

    // CRTP "inherited" functions
    void reserve(int, int);
    void reserveAdditional(int, int);
    int addParticle(const Particle&);
    void addParticle(const vector<Particle>&);
    int addParticle(const Particle&, Characteristic*);
    void removeAt(int);
    Particle makeParticle(int);
    vector<Particle> getParticles();
    void wrap(RealType&, RealType&);
    void wrap(RealType&);
    RealType getPhi();
    void updatePositionRecord();
    void setInitialPositions();

  protected:

    // Private CRTP "inherited" functions    
    inline void compressArrays();
    
  };

}

#endif // __SIM_DATA_HPP__
