#include "store_data.hpp"
// Other files
#include "../base/simdata.hpp"
#include "../base/topology.hpp"

namespace GFlowSimulation {

  StoreData::StoreData() : bounds(2), dataWidth(0), sim_dimensions(0), nTypes(0), write_processor_info(true) {}

  void StoreData::initialize(shared_ptr<SimData> sd) {
    // Make sure simData is non-null
    if (sd==nullptr) return;
    // Set simdata
    simData = sd;
    // Get basic data
    bounds = simData->getBounds();
    sim_dimensions = simData->getSimDimensions();
    nTypes = simData->ntypes();
    // Make sure there are types.
    if (nTypes==0) return;
    // Get data positions, calculate data width.
    dataWidth = 0;
    vector<string> temp;
    for (auto entry : vector_data_entries) {
      int pos = simData->getVectorData(entry);
      if (-1<pos) {
        temp.push_back(entry);
        vector_data_positions.push_back(pos);
        dataWidth += sim_dimensions;
      }
    }
    vector_data_entries = temp;
    temp.clear();
    for (auto entry : scalar_data_entries) {
      int pos = simData->getScalarData(entry);
      if (-1<pos) {
        temp.push_back(entry);
        scalar_data_positions.push_back(pos);
        ++dataWidth;
      }
    }
    scalar_data_entries = temp;
    temp.clear();
    for (auto entry : integer_data_entries) {
      int pos = simData->getIntegerData(entry);
      if (-1<pos) {
        temp.push_back(entry);
        integer_data_positions.push_back(pos);
        ++dataWidth;
      }
    }
    integer_data_entries = temp;

    // Processor data is takes up one integer slot.
    if (write_processor_info) ++dataWidth;
  }

  void StoreData::set_vector_data(const vector<string>& v) {
    vector_data_entries = v;
  }

  void StoreData::set_scalar_data(const vector<string>& v) {
    scalar_data_entries = v;
  }

  void StoreData::set_integer_data(const vector<string>& v) {
    integer_data_entries = v;
  }

  void StoreData::set_data_boundary(const Bounds& bnds) {
    bounds = bnds;
  }

  const vector<string>& StoreData::get_vector_data() const {
    return vector_data_entries;
  }

  const vector<string>& StoreData::get_scalar_data() const {
    return scalar_data_entries;
  }

  const vector<string>& StoreData::get_integer_data() const {
    return integer_data_entries;
  }

  void StoreData::store(vector<float>& data) {
    data.clear();
    int number = simData->number_owned();

    // If there are no particles/no data, return
    if (dataWidth==0 || simData==nullptr || number==0) return;
    // Set up data
    data = vector<float>(dataWidth*number, 0);

    // Fill the array
    int data_pointer = 0;
    for (int n=0; n<number; ++n) {
      // If not a particle, continue
      if (simData->Type(n)<0 || !bounds.contains(simData->X(n))) continue;
      // Copy data
      for (auto v : vector_data_positions) {
        auto vd = simData->VectorData(v);
        if (!vd.isnull()) copyVec(vd(n), &data[data_pointer], sim_dimensions);
        data_pointer += sim_dimensions;
      }
      for (auto s : scalar_data_positions) {
        auto sd = simData->ScalarData(s);
        if (!sd.isnull()) data[data_pointer] = sd[n];
        ++data_pointer;
      }
      for (auto i : integer_data_positions) {
        auto id = simData->IntegerData(i);
        if (!id.isnull()) data[data_pointer] = id[n]; 
        ++data_pointer;
      }
      if (write_processor_info) {
        if (simData->isReal(n)) data[data_pointer] = 2*simData->getTopology()->getRank();
        else data[data_pointer] = 2*simData->getTopology()->getRank() + 1;
        ++data_pointer;
      }
    }
  }

  bool StoreData::write(const string& fileName, const vector<vector<float> >& positions) {
    // Print data to csv
    ofstream fout(fileName);
    if (fout.fail()) return false;

    // Write the file header
    writeHeader(fout, positions.size());

    // Print out the actual data - first the number of particles, then the particle data
    for (auto &v : positions) fout << v.size() << "," << toCSV(v) << "\n";
    fout.close();

    // Return success
    return true;
  }

  bool StoreData::write(const string& fileName, const vector<float>& positions) {
    // Print data to csv
    ofstream fout(fileName);
    if (fout.fail()) return false;

    // Write the file header
    writeHeader(fout, 1);

    // Print out the actual data - first the number of particles, then the particle data
    fout << positions.size() << "," << toCSV(positions) << "\n";
    fout.close();

    // Return success
    return true;
  }

  const Bounds& StoreData::getBounds() const {
    return bounds;
  }

  int StoreData::getDataWidth() const {
    return dataWidth;
  }

  int StoreData::getDimensions() const {
    return sim_dimensions;
  }

  int StoreData::getNTypes() const {
    return nTypes;
  }

  inline void StoreData::writeHeader(ofstream& fout, int iters) const {
    // Print data width, dimensions, data iterations, ntypes
    fout << dataWidth << "," << sim_dimensions << "," << iters << "," << nTypes << "\n";

    // Print bounds - mins, then maxes
    for (int i=0; i<sim_dimensions; ++i) 
      fout << bounds.min[i] << ",";
    for (int i=0; i<sim_dimensions; ++i) {
      fout << bounds.max[i];
      if (i!=sim_dimensions-1) fout << ",";
    }
    fout << "\n";

    // Vector data types
    fout << vector_data_entries.size() << ",";
    for (int i=0; i<vector_data_entries.size(); ++i) {
      fout << vector_data_entries[i];
      if (i!=vector_data_entries.size()-1) fout << ",";
    }
    fout << "\n";

    // Scalar data types
    fout << scalar_data_entries.size() << ",";
    for (int i=0; i<scalar_data_entries.size(); ++i) {
      fout << scalar_data_entries[i];
      if (i!=scalar_data_entries.size()-1) fout << ",";
    }
    fout << "\n";
    
    // Integer data types
    int integer_data_size = integer_data_entries.size(); 
    if (write_processor_info) ++integer_data_size;
    fout << integer_data_size << ",";
    for (int i=0; i<integer_data_entries.size(); ++i) {
      fout << integer_data_entries[i];
      if (i!=integer_data_entries.size()-1) fout << ",";
    }
    // Processor info.
    if (write_processor_info) {
      if (integer_data_size>1) fout << ",Proc";
      else fout << "Proc";
    }
    fout << "\n";
    
  }

}
