#include "dataobjects/dataobjecttypes/data-base/particle-store-data.hpp"

using namespace GFlowSimulation;

void ParticleStoreData::add_vector_data_entry(const std::string& entry) {
  vector_data_entries.push_back(entry);
}

void ParticleStoreData::add_magnitude_data_entry(const std::string& entry) {
  magnitude_data_entries.push_back(entry);
}

void ParticleStoreData::add_scalar_data_entry(const std::string& entry) {
  scalar_data_entries.push_back(entry);
}

void ParticleStoreData::add_integer_data_entry(const std::string& entry) {
  integer_data_entries.push_back(entry);
}

void ParticleStoreData::clear_all_data_entries() {
  vector_data_entries.clear();
  scalar_data_entries.clear();
  integer_data_entries.clear();
}
