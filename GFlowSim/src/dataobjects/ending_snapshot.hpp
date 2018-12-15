#ifndef __ENDING_SNAPSHOT_HPP__GFLOW__
#define __ENDING_SNAPSHOT_HPP__GFLOW__

#include "../base/dataobject.hpp"

namespace GFlowSimulation {

  class EndingSnapshot : public DataObject {
  public:
    //! @brief Constructor
    EndingSnapshot(GFlow*);

    //! @brief Destructor.
    ~EndingSnapshot();

    //! @brief Store the initial positions of the particles.
    virtual void pre_integrate() override;

    virtual void post_integrate() override;

    //! @brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    // Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

  private:
    
    //! @brief Helper function that stores all the neccesary data in a vector.
    template<typename T> void store_data(vector<T>& data, vector<DataType>& d_types) {
      // Fill the array of data
      for (int i=0; i<Base::simData->size(); ++i) {
        if (Base::simData->Type(i)!=-1)
          for (const auto type : d_types) get_data(data, type, i);
      }
    }

    //! @brief Store a particular type of data.
    template<typename T> void get_data(vector<T>& data, DataType type, int id) {
    switch(type) {
      case DataType::POSITION: {
        for (int d=0; d<sim_dimensions; ++d) data.push_back(Base::simData->X(id, d));
        break;
      }
      case DataType::VELOCITY: {
        for (int d=0; d<sim_dimensions; ++d) data.push_back(Base::simData->V(id, d));
        break;
      }
      case DataType::SIGMA: {
        data.push_back(Base::simData->Sg(id));
        break;
      }
      case DataType::TYPE: {
        data.push_back(Base::simData->Type(id));
        break;
      }
      case DataType::DISTANCE: {
        gflow->getDisplacement(Base::simData->X(id), &initial_data.at(id*sim_dimensions), vdata);
        data.push_back( magnitudeVec(vdata, sim_dimensions) );
        break;
      }
    }
  }

    //! @brief Contains all the relevant data for a time step.
    //!
    //! Each time step contains a length ( [number of particles] * [dataWidth] ) array 
    //! for the positions of the particles in [DIMENSIONS] dimensions, and any other data we want to store.
    vector<double> final_data; 

    //! @brief Initial positions of particle.
    vector<RealType> initial_data;

    //! @brief The amount of data we collect per particle.
    int dataWidth;

    //! @brief A holder for vector data
    RealType *vdata;

    vector<DataType> data_types;
  };

}
#endif // __ENDING_SNAPSHOT_HPP__GFLOW__