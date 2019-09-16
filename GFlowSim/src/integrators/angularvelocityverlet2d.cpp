#include "angularvelocityverlet2d.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../utility/vectormath.hpp"

namespace GFlowSimulation {

  AngularVelocityVerlet2d::AngularVelocityVerlet2d(GFlow *gflow) : Integrator(gflow) {};

  void AngularVelocityVerlet2d::pre_forces() {
    // Start the timer
    timer.start();

    // Call base class
    Integrator::pre_forces();

    // Look for data entries.
    RealType *sg = simData->Sg();
    RealType *im = simData->Im();
    int th_add = simData->getScalarData("Th");
    int om_add = simData->getScalarData("Om");
    int tq_add = simData->getScalarData("Tq");
    RealType *th = simData->ScalarData(th_add);
    RealType *om = simData->ScalarData(om_add);
    RealType *tq = simData->ScalarData(tq_add);
    // Check to make sure we have all the data we need.
    if (sim_dimensions!=2 || sg==nullptr || im==nullptr || om==nullptr || tq==nullptr) return;
    
    // --- First half kick

    // Number of (real - non ghost) particles
    const int total = simData->size_owned();
    if (total==0) return;
    // Half a timestep
    RealType hdt = 0.5*gflow->getDT();
    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<total; i+=2) {
      RealType iII = 2*im[i]/sqr(sg[i]); // Inverse moment of inertia.
      om[i] += hdt*iII*tq[i];
    }
    #else
    // Put hdt into a simd vector
    simd_float _2hdt = simd_set1(2*hdt);
    int i;
    for (i=0; i<total-simd_data_size; i+=simd_data_size) {

      simd_float torque = simd_load(&tq[i]);
      simd_float omega  = simd_load(&om[i]);
      simd_float invM   = simd_load(&im[i]);
      simd_float sigma  = simd_load(&sg[i]);

      simd_float hII = invM/(sigma*sigma);
      simd_float vec = _2hdt*hII*torque;
      omega += vec;

      simd_store(vec, &om[i]);
    }
    // Left overs
    for (; i<total; ++i) {
      RealType II = 2*im[i]/sqr(sg[i]); // Inverse moment of inertia.
      om[i] += hdt*II*tq[i];
    }
    #endif

    if (th==nullptr) {
      timer.stop();
      return;
    }

    // --- Update angles
    // Do serially
    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<total; ++i) {
      th[i] += dt*om[i];
    }
    // Do with SIMD
    #else
    // Set dt
    simd_float dt_vec = simd_set1(dt);
    for (i=0; i<=total-simd_data_size; i+=simd_data_size) {
      simd_float Th = simd_load(&th[i]);
      simd_float Om = simd_load(&om[i]);
      simd_float dTh = Om*dt_vec;
      simd_float Th_new = Th+dTh;
      simd_store(Th_new, &th[i]);
    }
    // Left overs
    for (; i<total; ++i) {
      th[i] += dt*om[i];
    }
    #endif

    // Stop timer
    timer.stop();
  }

  void AngularVelocityVerlet2d::post_forces() {
    // Start the timer
    timer.start();

    // Call to parent class
    Integrator::post_forces();

    // Look for data entries.
    RealType *sg = simData->Sg();
    RealType *im = simData->Im();
    int om_add = simData->getScalarData("Om");
    int tq_add = simData->getScalarData("Tq");
    RealType *om = simData->ScalarData(om_add);
    RealType *tq = simData->ScalarData(tq_add);
    // Check to make sure we have all the data we need.
    if (sim_dimensions!=2 || sg==nullptr || im==nullptr || om==nullptr || tq==nullptr) return;
    
    // --- Second half kick

    // Number of (real - non ghost) particles
    const int total = simData->size_owned();
    if (total==0) return;
    // Half a timestep
    RealType hdt = 0.5*gflow->getDT();

    #if SIMD_TYPE==SIMD_NONE
    for (int i=0; i<total; i+=2) {
      RealType II = 2*im[i]/sqr(sg[i]); // Inverse moment of inertia.
      om[i] += hdt*II*tq[i];
    }
    #else
    // Put hdt into a simd vector
    simd_float _2hdt = simd_set1(2*hdt);
    int i;
    for (i=0; i<total-simd_data_size; i+=simd_data_size) {

      simd_float torque = simd_load(&tq[i]);
      simd_float omega  = simd_load(&om[i]);
      simd_float invM   = simd_load(&im[i]);
      simd_float sigma  = simd_load(&sg[i]);

      simd_float hII = invM/(sigma*sigma);
      simd_float vec = _2hdt*hII*torque;
      omega += vec;

      simd_store(vec, &om[i]);
    }
    // Left overs
    for (; i<total; ++i) {
      RealType II = 2*im[i]/sqr(sg[i]); // Inverse moment of inertia.
      om[i] += hdt*II*tq[i];
    }
    #endif

    // Stop timer
    timer.stop();
  }

}