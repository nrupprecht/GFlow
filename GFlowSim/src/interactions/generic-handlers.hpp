#ifndef __GENERIC_HANDLERS_HPP__GFLOW__
#define __GENERIC_HANDLERS_HPP__GFLOW__

#include "../base/interaction.hpp"
#include "../utility/generic-dimension.hpp"

namespace GFlowSimulation {

  template<int dims, class ForceType> class VerletList : public Interaction {
  public:
    //! \brief Constructor.
    VerletList(GFlow *gflow) : Interaction(gflow) {};

    //! \brief Add a pair of particles whose distance may need to be wrapped.
    //!
    //! Adds the pair to the verlet_wrap list.
    virtual void addPair(const int id1, const int id2) override {
      if (id1!=last_head_wrap) {
        first_neighbor_wrap.push_back(verlet_wrap.size());
        head_list_wrap.push_back(id1);
        last_head_wrap = id1;
      }
      verlet_wrap.push_back(id2);
      // Increment pairs.
      ++pairs;
    }
    
    //! \brief Add a pair of interacting particles that do not need to have their distances wrapped.
    //!
    //! Adds the pair to the verlet_wrap list.
    virtual void addPairNW(const int id1, const int id2) override {
      if (id1!=last_head) {
        first_neighbor.push_back(verlet.size());
        head_list.push_back(id1);
        last_head = id1;
      }
      verlet.push_back(id2);
      // Increment pairs.
      ++pairs;
    }

    virtual void close() override {
      // Mark the end of the last neighbor lists.
      first_neighbor_wrap.push_back(verlet_wrap.size());
      first_neighbor.push_back(verlet.size());
      // Reset heads.
      last_head = last_head_wrap = -1;
    }

    //! \brief Clears the verlet list.
    virtual void clear() override {
      // Clear verlet lists.
      Interaction::clear();
      // Clear first neighbor lists.
      first_neighbor.clear();
      first_neighbor_wrap.clear();
      // Clear head lists.
      head_list.clear();
      head_list_wrap.clear();
      // Reset other data.
      last_head = last_head_wrap = -1;
      pairs = 0;
    }

    //! \brief Returns the size of the verlet list.
    virtual int size() const override {
      return pairs;
    }

    //! \brief Iterate through all 
    virtual void interact() const override {
      // Common tasks 
      Interaction::interact();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Get the data pointers.
      RealType **x = simData->X();
      RealType **f = simData->F();
      RealType *rd = simData->Sg();
      int    *type = simData->Type();

      // Make sure all needed pointers are non null.
      // \todo Should probably have some sort of global error message system.
      if (x==nullptr || f==nullptr || rd==nullptr || type==nullptr) return;

      // Needed constants
      RealType R1, R2, rsqr, r, X[dims];

      // --- Go through all particles in verlet.
      for (int i=0; i<head_list.size(); ++i) {
        // Get id of the head particle of this list.
        int id1 = head_list[i];
        // Check if the type of the head is good.
        if (type[id1]<0) continue;
        // Iterate through head particle's neighbors.
        for (int j=first_neighbor[i]; j<first_neighbor[i+1]; ++j) {
          // Get id of neighbor.
          int id2 = verlet[j];
          // Check if the type of the neighbor is good.
          if (type[id2]<0) continue;
          // Calculate displacement.
          subtract_vec<dims>(x[id1], x[id2], X);
          // Calculate squared distance
          rsqr = dot_vec<dims>(X, X);
          // Get radii
          R1 = rd[id1];
          R2 = rd[id2];
          // If close, interact.
          if (rsqr < sqr((R1 + R2)*cutoff)) 
            static_cast<const ForceType*>(this)->force(id1, id2, R1, R2, rsqr, X, f);
        }
      }

      // --- Go through all particles in verlet wrap.
      if (verlet_wrap.empty()) return;

      // Get the bounds and boundary conditions
      vector<BCFlag> bcs;
      vector<RealType> widths;
      // Fill the vectors.
      for (int d=0; d<dims; ++d) {
        bcs.push_back(gflow->getBC(d));
        widths.push_back(gflow->getBounds().wd(d));
      }

      // --- Go through all particles in verlet.
      for (int i=0; i<head_list_wrap.size(); ++i) {
        // Get id of the head particle of this list.
        int id1 = head_list_wrap[i];
        // Iterate through head particle's neighbors.
        for (int j=first_neighbor_wrap[i]; j<first_neighbor_wrap[i+1]; ++j) {
          // Get id of neighbor.
          int id2 = verlet_wrap[j];

          // Check if the types are good
          if (type[id1]<0 || type[id2]<0) continue;
          // Calculate displacement.
          subtract_vec<dims>(x[id1], x[id2], X);
          // Harmonic corrections to distance.
          harmonic_correction<dims>(bcs.data(), X, widths.data());
          // Calculate squared distance
          rsqr = dot_vec<dims>(X, X);
          // Get radii
          R1 = rd[id1];
          R2 = rd[id2];
          // If close, interact.
          if (rsqr < sqr((R1 + R2)*cutoff)) 
            static_cast<const ForceType*>(this)->force(id1, id2, R1, R2, rsqr, X, f);
        }
      }
    }

  private:

    int last_head_wrap = -1;
    int last_head = -1;

    vector<int> first_neighbor;
    vector<int> first_neighbor_wrap;

    vector<int> head_list;
    vector<int> head_list_wrap;

    //! \brief The number of interaction pairs.
    int pairs = 0;
  };


  
  template<int dims, class ForceType> class VerletListVecPairs : public Interaction {
  public:
    //! \brief Constructor.
    VerletListVecPairs(GFlow *gflow) : Interaction(gflow) {};

    //! \brief Iterate through all 
    virtual void interact() const override {
      // Common tasks 
      Interaction::interact();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Get the data pointers.
      RealType **x = simData->X();
      RealType **f = simData->F();
      RealType *rd = simData->Sg();
      int    *type = simData->Type();

      // Make sure all needed pointers are non null.
      // \todo Should probably have some sort of global error message system.
      if (x==nullptr || f==nullptr || rd==nullptr || type==nullptr) return;

      // Needed constants
      RealType R1, R2, rsqr, r, X[dims];

      // --- Go through all particles in verlet.
      for (int i=0; i<verlet.size(); i+=2) {
        int id1 = verlet[i];
        int id2 = verlet[i+1];
        // Check if the types are good
        if (type[id1]<0 || type[id2]<0) continue;
        // Calculate displacement.
        subtract_vec<dims>(x[id1], x[id2], X);
        // Calculate squared distance
        rsqr = dot_vec<dims>(X, X);
        // Get radii
        R1 = rd[id1];
        R2 = rd[id2];
        // If close, interact.
        if (rsqr < sqr((R1 + R2)*cutoff)) 
          static_cast<const ForceType*>(this)->force(id1, id2, R1, R2, rsqr, X, f);
      }

      // --- Go through all particles in verlet wrap.
      if (verlet_wrap.empty()) return;

      // Get the bounds and boundary conditions
      vector<BCFlag> bcs;
      vector<RealType> widths;
      // Fill the vectors.
      for (int d=0; d<dims; ++d) {
        bcs.push_back(gflow->getBC(d));
        widths.push_back(gflow->getBounds().wd(d));
      }

      // --- Go through all particles in verlet_wrap.
      for (int i=0; i<verlet_wrap.size(); i+=2) {
        int id1 = verlet_wrap[i];
        int id2 = verlet_wrap[i+1];
        // Check if the types are good
        if (type[id1]<0 || type[id2]<0) continue;
        // Calculate displacement.
        subtract_vec<dims>(x[id1], x[id2], X);
        // Harmonic corrections to distance.
        harmonic_correction<dims>(bcs.data(), X, widths.data());
        // Calculate squared distance
        rsqr = dot_vec<dims>(X, X);
        // Get radii
        R1 = rd[id1];
        R2 = rd[id2];
        // If close, interact.
        if (rsqr < sqr((R1 + R2)*cutoff)) 
          static_cast<const ForceType*>(this)->force(id1, id2, R1, R2, rsqr, X, f);
      }
    }
  };



  template<int dims, class ForceType> class VerletListPairs : public Interaction {
  public:
    //! \brief Constructor.
    VerletListPairs(GFlow *gflow) : Interaction(gflow) {};

    //! \brief Add a pair of particles whose distance may need to be wrapped.
    //!
    //! Adds the pair to the verlet_wrap list.
    virtual void addPair(const int id1, const int id2) override {
      verlet_list_wrap.push_back(std::make_pair(id1, id2));
    }
    
    //! \brief Add a pair of interacting particles that do not need to have their distances wrapped.
    //!
    //! Adds the pair to the verlet_wrap list.
    virtual void addPairNW(const int id1, const int id2) override {
      verlet_list.push_back(std::make_pair(id1, id2));
    }

    //! \brief Clears the verlet list.
    virtual void clear() override {
      verlet_list.clear();
      verlet_list_wrap.clear();
    }

    //! \brief Returns the size of the verlet list.
    virtual int size() const override {
      return verlet_list.size() + verlet_list_wrap.size();
    }

    //! \brief Iterate through all 
    virtual void interact() const override {
      // Common tasks 
      Interaction::interact();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Get the data pointers.
      RealType **x = simData->X();
      RealType **f = simData->F();
      RealType *rd = simData->Sg();
      int    *type = simData->Type();

      // Make sure all needed pointers are non null.
      // \todo Should probably have some sort of global error message system.
      if (x==nullptr || f==nullptr || rd==nullptr || type==nullptr) return;

      // Needed constants
      RealType R1, R2, rsqr, r, X[dims];

      // --- Go through all particles in verlet.
      for (auto vpair : verlet_list) {
        int id1 = vpair.first;
        int id2 = vpair.second;
        // Check if the types are good
        if (type[id1]<0 || type[id2]<0) continue;
        // Calculate displacement.
        subtract_vec<dims>(x[id1], x[id2], X);
        // Calculate squared distance
        rsqr = dot_vec<dims>(X, X);
        // Get radii
        R1 = rd[id1];
        R2 = rd[id2];
        // If close, interact.
        if (rsqr < sqr((R1 + R2)*cutoff)) 
          static_cast<const ForceType*>(this)->force(id1, id2, R1, R2, rsqr, X, f);
      }

      // --- Go through all particles in verlet wrap.
      if (verlet_list_wrap.empty()) return;

      // Get the bounds and boundary conditions
      vector<BCFlag> bcs;
      vector<RealType> widths;
      // Fill the vectors.
      for (int d=0; d<dims; ++d) {
        bcs.push_back(gflow->getBC(d));
        widths.push_back(gflow->getBounds().wd(d));
      }

      // --- Go through all particles in verlet_wrap.
      for (auto vpair : verlet_list_wrap) {
        int id1 = vpair.first;
        int id2 = vpair.second;
        // Check if the types are good
        if (type[id1]<0 || type[id2]<0) continue;
        // Calculate displacement.
        subtract_vec<dims>(x[id1], x[id2], X);
        // Harmonic corrections to distance.
        harmonic_correction<dims>(bcs.data(), X, widths.data());
        // Calculate squared distance
        rsqr = dot_vec<dims>(X, X);
        // Get radii
        R1 = rd[id1];
        R2 = rd[id2];
        // If close, interact.
        if (rsqr < sqr((R1 + R2)*cutoff)) 
          static_cast<const ForceType*>(this)->force(id1, id2, R1, R2, rsqr, X, f);
      }
    }

  private:
    //! \brief The no-wrap verlet list. A list of pairs of interacting particles.
    vector<pair<int, int> > verlet_list;

    //! \brief The wrap verlet list. A list of pairs of interacting particles.
    vector<pair<int, int> > verlet_list_wrap;
  };

}
#endif // __GENERIC_HANDLERS_HPP__GFLOW__