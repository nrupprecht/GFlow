#include "field_properties.hpp"
// Other files
#include "load_data.hpp"
#include "../utility/vectormath.hpp"
#include <fstream>

namespace GFlowSimulation {

  void FieldProperties::create_field(const vector<double>& data, const Bounds& bounds, int entries, int dat_width) {
    // Set data width
    data_width = dat_width;
    // Assume 2 dimensions
    setPlaces(2);

    // Resolution
    double res = 0.15;

    int bx = static_cast<int>(bounds.wd(0)/res);
    int by = static_cast<int>(bounds.wd(1)/res);

    // Check for good bounds
    if (bx<1 || by<1) return;

    // Create a binning
    Binning2d binning(bx, by);
    binning.setBounds(bounds);
    binning.bin_data(data, entries, data_width);

    // Set up fields
    phi_field.setSize(bx, by);
    phi_field.setBounds(bounds);
    ave_v_field.setSize(bx, by);
    ave_v_field.setBounds(bounds);
    var_v_field.setSize(bx, by);
    var_v_field.setBounds(bounds);
    ave_v_vector_field.setSize(bx, by);
    ave_v_vector_field.setBounds(bounds);

    // Stencil
    int sx = 1, sy = 1;
    double ave_vx = 0, ave_vy = 0, vsqr = 0, var = 0;
    int num = 0;
    for (int y=0; y<by; ++y) {
      for (int x=0; x<bx; ++x) {
        // Reset values
        ave_vx = ave_vy = vsqr = var = 0;
        num = 0;
        // Calculate average velocity
        num = binning.at(x, y).size();
        for (const auto& id : binning.at(x, y)) {
          double vx = vel(id, 0, data);
          double vy = vel(id, 1, data);
          vsqr += (sqr(vx) + sqr(vy));
          ave_vx += vx; ave_vy += vy;
        }
        // Divide, to get average
        ave_vx = num>0 ? ave_vx/num : 0;
        ave_vy = num>0 ? ave_vy/num : 0;
        // Find variance
        for (const auto& id : binning.at(x, y))
          var += (sqr(vel(id, 0, data) - ave_vx) + sqr(vel(id, 1, data) - ave_vy));

        // Update fields
        phi_field.at(x,y)   = num;
        ave_v_field.at(x,y) = num>0 ? sqrt(vsqr)/num : 0;
        var_v_field.at(x,y) = num>1 ? sqrt(var)/(num-1) : 0;
        ave_v_vector_field.at(x,y) = vec2d(ave_vx, ave_vy);
      }
    }
  }

  void FieldProperties::load_and_create(string directory) {
    LoadData loader;
    // Load data
    if (loader.load(directory+"/data.csv")) {
      const auto& data = loader.getData();
      for (int i=0; i<data.size(); ++i) {
        // Create the field of data
        create_field(data[i], loader.getBounds(), loader.getEntries(i), loader.getDataWidth());
        // Save the field of data
        std::ofstream fout(directory+"/aveV-VF-"+toStr(i)+".csv");
        if (!fout.fail()) 
          fout << toCSV_aveVVectorField();
        fout.close();
      }
    }
  }

  string FieldProperties::toStrMMA_phiField() {
    return phi_field.toMMA();
  }

  string FieldProperties::toStrMMA_aveVField() {
    return ave_v_field.toMMA();
  }

  string FieldProperties::toStrMMA_varVField() {
    return var_v_field.toMMA();
  }

  string FieldProperties::toCSV_phiField() {
    return phi_field.toCSV();
  }

  string FieldProperties::toCSV_aveVField() {
    return ave_v_field.toCSV();
  }

  string FieldProperties::toCSV_varVField() {
    return var_v_field.toCSV();
  }

  string FieldProperties::toCSV_aveVVectorField() {
    return ave_v_vector_field.toCSV();
  }

  inline double FieldProperties::vel(int id, int d, const vector<double>& data) {
    return data.at(data_width*id + vel_place + d);
  }

  inline void FieldProperties::setPlaces(const int d) {
    pos_place = 0;
    vel_place = d;
    sg_place = 2*d;
    type_place = 2*d+1; 
    distance_place = 2*d+2;
  }

}