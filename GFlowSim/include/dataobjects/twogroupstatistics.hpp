#ifndef __TWO_GROUP_STATISTICS_HPP__GFLOW__
#define __TWO_GROUP_STATISTICS_HPP__GFLOW__

#include "../base/dataobject.hpp"
#include "../other/group.hpp"

namespace GFlowSimulation {

  class TwoGroupStatistics : public DataObject {
  public:
    //! \brief Default constructor.
    TwoGroupStatistics(GFlow*);

    //! \brief Clear out any old data.
    virtual void pre_integrate() override;

    //! \brief Collect the position data from simdata --- happens during the post-step phase
    virtual void post_step() override;

    //! \brief Write data to a file - if true, the string is a path, and you should use your own name as the file name.
    // Returns true for success.
    virtual bool writeToFile(string, bool=true) override;

    //! \brief Set the number of x bins.
    void setBins(int, int);

    //! \brief Give two groups to be the groups.
    void setGroups(const Group&, const Group&);

  private:
    //! \brief The two groups of atoms
    Group group1, group2;

    //! \brief The histogram of distance, relative velocity.
    vector<vector<int> > count_v;
    //! \brief The histogram of distance, relative force.
    vector<vector<int> > count_f;

    //! \brief The number of x (distance) bins to use.
    int bins_x;

    //! \brief The number of y (velocity or force) bins to use.
    int bins_y;

    //! \brief The cutoff distance.
    RealType r_max = 5.;

    //! \brief The max absolute value of velocity to bin.
    RealType v_max = 0.25;

    //! \brief The max absolute value of force to bin.
    RealType f_max = 0.25;
  };

}
#endif // __TWO_GROUP_STATISTICS_HPP__GFLOW__