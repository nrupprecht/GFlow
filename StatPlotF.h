#ifndef __STAT_PLOT_F_H__
#define __STAT_PLOT_F_H__

typedef void (*StatPlotF) (const vector<PData>&, vector<vec2>&, const double, const double);

inline void Plot_Force_Vs_Depth(const vector<PData>& positions, vector<vec2>& statVector, const double lower, const double upper) {
  int bins = statVector.size();
  if (bins==0) return;
  // Bin data
  vector<double> pressBins(bins, 0);
  vector<int> counts(bins, 0);
  double dq = (upper-lower)/bins;
  for (const auto &p : positions) {
    vec2 pos = std::get<0>(p);
    int b = (upper-pos.y)/dq;
    pressBins.at(b) += std::get<4>(p);
    ++counts.at(b);
  }
  for (int i=0; i<bins; ++i) {
    statVector.at(i).y += counts.at(i)>0 ? pressBins.at(i)/static_cast<double>(counts.at(i)) : 0;
  }
}

inline void Plot_Force_Sigma_Vs_Depth(const vector<PData>& positions, vector<vec2>& statVector, const double lower, const double upper) {
  int bins = statVector.size();
  if (bins==0) return;
  // Bin data
  vector<vector<double> > forces(bins);
  double dq = (upper-lower)/bins;
  for (const auto &p : positions) {
    vec2 pos = std::get<0>(p);
    int b = (upper-pos.y)/dq;
    forces.at(b).push_back(std::get<4>(p));
  }
  for (int i=0; i<bins; ++i) {
    int num = forces.at(i).size();
    if (num>0) {
      // Find average
      double ave = 0;
      for (const auto &f : forces.at(i)) ave += f;
      ave /= num;
      // Find standard deviation
      double std = 0;
      for (const auto &f : forces.at(i)) std += sqr(ave-f);
      std /= num;
      statVector.at(i).y += std;
    }
    else; // Add zero = do nothing
  }
}

#endif // __STAT_PLOT_F_H__
