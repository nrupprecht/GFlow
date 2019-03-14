#ifndef __PARTICLE_CONTAINER_HPP__GFLOW__
#define __PARTICLE_CONTAINER_HPP__GFLOW__

namespace GFlowSimulation {

  class ParticleContainer : public Base {
  public:

    ParticleContainer(GFlow*);

    virtual void initialize();



  private:

    RealType *data;

    std::map<string, int> vector_data_offsets;
    std::map<string, int> scalar_data_offsets;
    std::map<string, int> integer_data_offsets;

  };

}
#endif // __PARTICLE_CONTAINER_HPP__GFLOW__