#ifndef __GENERIC_HANDLERS_HPP__GFLOW__
#define __GENERIC_HANDLERS_HPP__GFLOW__

#include "../base/interaction.hpp"
#include "../utility/generic-dimension.hpp"

#include "../base/topology.hpp"

namespace GFlowSimulation {

  template<int dims, class ForceType, bool using_history> class VerletList : public Interaction {
  public:
    //! \brief Constructor.
    VerletList(GFlow *gflow) : Interaction(gflow) {
      static_assert(using_history==false, "VerletList cannot work with interaction history.");
      last_head[0] = last_head[1] = last_head[2] = -1;
    };

    //! \brief Add a pair of particles whose distance may need to be wrapped.
    //!
    //! Adds the pair to the verlet_wrap list.
    virtual void addPair(const int id1, const int id2, const int list) override {
      if (id1!=last_head[list]) {
        first_neighbor[list].push_back(verlet[list].size());
        head_list[list].push_back(id1);
        last_head[list] = id1;
      }
      verlet[list].push_back(id2);
      // Increment pairs.
      ++pairs;
    }

    virtual void close() override {
      // Mark the end of the last neighbor lists.
      first_neighbor[0].push_back(verlet[0].size());
      first_neighbor[1].push_back(verlet[1].size());
      first_neighbor[2].push_back(verlet[2].size());
      // Reset heads.
      last_head[0] = last_head[1] = last_head[2] = -1;
    }

    //! \brief Clears the verlet list.
    virtual void clear() override {
      // Clear verlet lists.
      Interaction::clear();
      // Clear first neighbor lists.
      first_neighbor[0].clear();
      first_neighbor[1].clear();
      first_neighbor[2].clear();
      // Clear head lists.
      head_list[0].clear();
      head_list[1].clear();
      head_list[2].clear();
      // Reset heads.
      last_head[0] = last_head[1] = last_head[2] = -1;
      // Reset number of pairs.
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

      // Do forces.
      dolist<false>(0);
      dolist<true>(1);
    }

    //! \brief Iterate through all 
    virtual void interact_ghosts() const override {
      // Common tasks 
      Interaction::interact_ghosts();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Do forces for ghosts.
      potential_factor = 0.5f;
      dolist<true, 1>(2);
      potential_factor = 1.0f;
    }

  private:

    //! \brief Do the forces for a single list. 
    template<bool wrapping, int first_particle_type=0> inline void dolist(const int list) const {
      // Reference to the requested list.
      const vector<int> &verlet_list = verlet[list];
      // If there are no particles, return.
      if (verlet_list.empty()) return;
      // Get references to the rest of the list data we need.
      auto& first = first_neighbor[list];
      auto& heads = head_list[list];

      // Get the data pointers.
      auto x0 = simData->X<first_particle_type>();
      auto x1 = simData->X();
      auto rd0 = simData->Sg<first_particle_type>();
      auto rd1 = simData->Sg();

      // Set up internal arrays.
      static_cast<const ForceType*>(this)->template set_types<first_particle_type, 0>();

      // Needed constants
      RealType R1, R2, rsqr, X[dims];

      // Get the bounds and boundary conditions...
      BCFlag bcs[dims];
      RealType widths[dims];
      // real min_wrappable_distance;
      // ... but only if wrapping.
      if (wrapping) {
        // min_wrappable_distance = 0.5f*gflow->getBounds().wd(0);
        // Fill the vectors.
        for (int d=0; d<dims; ++d) {
          bcs[d] = gflow->getBC(d);
          widths[d] = gflow->getBounds().wd(d);
          // if (0.5f*widths[d]<min_wrappable_distance) min_wrappable_distance = 0.5f*widths[d];
        }
      }

      // --- Go through all particles in verlet.
      for (int i=0; i<heads.size(); ++i) {
        // Get id of the head particle of this list.
        int id1 = heads[i];
        // Iterate through head particle's neighbors.
        for (int j=first[i]; j<first[i+1]; ++j) {
          // Get id of neighbor.
          int id2 = verlet_list[j];
          // Calculate displacement.
          subtract_vec<dims>(x0(id1), x1(id2), X);
          // Harmonic corrections to distance.
          if (wrapping) harmonic_correction<dims>(bcs, X, widths);
          // Calculate squared distance
          rsqr = dot_vec<dims>(X, X);
          // Get radii
          R1 = rd0(id1);
          R2 = rd1(id2);
          // If close, interact.
          if (rsqr < sqr((R1 + R2)*cutoff)) 
            static_cast<const ForceType*>(this)->force(id1, id2, R1, R2, rsqr, X);
        }
      }
    }

    //! \brief The most recent (and current) head particle for the list that we are adding particles to.
    //!
    //! If a pair of particles (id1, id2) is inserted into list n where id1!=last_head[n], then we start
    //! a new list, and update last_head[n] = id1.
    int last_head[3];

    //! \brief The (local) particle id of the i-th head.
    vector<int> head_list[3];

    //! \brief The first neighbor of the i-th head particle (for the n-th list) is located at first_neighbor[n][i]
    vector<int> first_neighbor[3];

    //! \brief The number of interaction pairs (across all the verlet lists in this interaction).
    int pairs = 0;
  };


  
  template<int dims, class ForceType, bool using_history> class VerletListVecPairs : public Interaction {
  public:
    //! \brief Constructor.
    VerletListVecPairs(GFlow *gflow) : Interaction(gflow) {
      static_assert(using_history==false, "VerletListVecPairs cannot work with interaction history.");
    };

    //! \brief Iterate through all 
    virtual void interact() const override {
      // Common tasks 
      Interaction::interact();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Do forces.
      dolist<false>(verlet[0]);
      dolist<true>(verlet[1]);
    }

    //! \brief Iterate through all 
    virtual void interact_ghosts() const override {
      // Common tasks 
      Interaction::interact_ghosts();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Do forces for ghosts.
      potential_factor = 0.5f;
      dolist<true, 1>(verlet[2]);
      potential_factor = 1.0f;
    }

  private:

    //! \brief Do the forces for a single list. 
    template<bool wrapping, int first_particle_type=0> inline void dolist(const vector<int>& verlet_list) const {
      // If the list is empty, return.
      if (verlet_list.empty()) return;

      // Get the data pointers.
      auto x0 = simData->X<first_particle_type>();
      auto x1 = simData->X();
      auto rd0 = simData->Sg<first_particle_type>();
      auto rd1 = simData->Sg();

      // Set up internal arrays.
      static_cast<const ForceType*>(this)->template set_types<first_particle_type, 0>();

      // Needed constants
      RealType R1, R2, rsqr, X[dims];

      // Get the bounds and boundary conditions...
      BCFlag bcs[dims];
      RealType widths[dims];
      // ... but only if wrapping.
      if (wrapping) {
        // Fill the vectors.
        for (int d=0; d<dims; ++d) {
          bcs[d] = gflow->getBC(d);
          widths[d] = gflow->getBounds().wd(d);
        }
      }

      // --- Go through all particles in verlet_list.
      for (int i=0; i<verlet_list.size(); i+=2) {
        int id1 = verlet_list[i];
        int id2 = verlet_list[i+1];
        // Calculate displacement.
        subtract_vec<dims>(x0(id1), x1(id2), X);
        // Harmonic corrections to distance.
        if (wrapping) harmonic_correction<dims>(bcs, X, widths);
        // Calculate squared distance
        rsqr = dot_vec<dims>(X, X);
        // Get radii
        R1 = rd0(id1);
        R2 = rd1(id2);
        // If close, interact.
        if (rsqr < sqr((R1 + R2)*cutoff)) 
          static_cast<const ForceType*>(this)->force(id1, id2, R1, R2, rsqr, X);
      }
    }
  };



  template<int dims, class ForceType, bool using_history> class VerletListPairs : public Interaction {
  public:
    //! \brief Constructor.
    VerletListPairs(GFlow *gflow) : Interaction(gflow) {
      static_assert(using_history==false, "VerletListPairs cannot work with interaction history.");
    };

    //! \brief Add a pair of particles whose distance may need to be wrapped.
    //!
    //! Adds the pair to the verlet_wrap list.
    virtual void addPair(const int id1, const int id2, const int list) override {
      verlet_pairs[list].push_back(std::make_pair(id1, id2));
    }

    //! \brief Clears the verlet list.
    virtual void clear() override {
      verlet_pairs[0].clear();
      verlet_pairs[1].clear();
      verlet_pairs[2].clear();
    }

    //! \brief Returns the size of the verlet list.
    virtual int size() const override {
      return verlet_pairs[0].size() + verlet_pairs[1].size() + verlet_pairs[2].size();
    }

    //! \brief Iterate through all 
    virtual void interact() const override {
      // Common tasks 
      Interaction::interact();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Do forces.
      dolist<false>(verlet_pairs[0]);
      dolist<true>(verlet_pairs[1]);
    }

    //! \brief Iterate through all 
    virtual void interact_ghosts() const override {
      // Common tasks 
      Interaction::interact_ghosts();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Do forces for ghosts.
      potential_factor = 0.5f;
      dolist<true, 1>(verlet_pairs[2]);
      potential_factor = 1.0f;
    }

  private:

    //! \brief Do the forces for a single list. 
    template<bool wrapping, int first_particle_type=0> inline void dolist(const vector<pair<int, int> >& verlet_list) const {
      // If the list is empty, return.
      if (verlet_list.empty()) return;

      // Get the data pointers.
      auto x0 = simData->X<first_particle_type>();
      auto x1 = simData->X();
      auto rd0 = simData->Sg<first_particle_type>();
      auto rd1 = simData->Sg();

      // Set up internal arrays.
      static_cast<const ForceType*>(this)->template set_types<first_particle_type, 0>();

      // Needed constants
      RealType R1, R2, rsqr, X[dims];

      // Get the bounds and boundary conditions...
      BCFlag bcs[dims];
      RealType widths[dims];
      // ... but only if wrapping.
      if (wrapping) {
        // Fill the vectors.
        for (int d=0; d<dims; ++d) {
          bcs[d] = gflow->getBC(d);
          widths[d] = gflow->getBounds().wd(d);
        }
      }

      // --- Go through all particles in verlet_wrap.
      for (const auto vpair : verlet_list) {
        int id1 = vpair.first;
        int id2 = vpair.second;
        // Calculate displacement.
        subtract_vec<dims>(x0(id1), x1(id2), X);
        // Harmonic corrections to distance.
        if (wrapping) harmonic_correction<dims>(bcs, X, widths);
        // Calculate squared distance
        rsqr = dot_vec<dims>(X, X);
        // Get radii
        R1 = rd0(id1);
        R2 = rd1(id2);
        // If close, interact.
        if (rsqr < sqr((R1 + R2)*cutoff))
          static_cast<const ForceType*>(this)->force(id1, id2, R1, R2, rsqr, X);
      }
    }

    //! \brief The verlet lists.
    vector<pair<int, int> > verlet_pairs[3];
  };


  template<int dims, class ForceType, bool using_history> class VerletListVector: public Interaction {
  public:
    //! \brief Constructor.
    VerletListVector(GFlow *gflow) : Interaction(gflow) {};

    //! \brief Add a pair of particles whose distance may need to be wrapped.
    //!
    //! Adds the pair to the verlet_wrap list.
    virtual void addPair(const int id1, const int id2, const int list) override {
      id1<id2 ? neighbor_list[list][id1].push_back(id2) : neighbor_list[list][id2].push_back(id1);
      ++_size;
    }

    //! \brief Clears the verlet list.
    virtual void clear() override {
      for (int i=0; i<3; ++i)
        for (auto& vec : neighbor_list[i]) vec.clear();
      // Reset size
      _size = 0;
    }

    //! \brief Returns the size of the verlet list.
    virtual int size() const override {
      return _size;
    }

    //! \brief Iterate through all 
    virtual void interact() const override {
      // Common tasks 
      Interaction::interact();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Do forces.
      dolist<false>(neighbor_list[0]);
      dolist<true>(neighbor_list[1]);
    }

    //! \brief Iterate through all 
    virtual void interact_ghosts() const override {
      // Common tasks 
      Interaction::interact_ghosts();

      // Do dimensional check.
      // \todo Should probably have some sort of global error message system.
      if (sim_dimensions!=dims) return;

      // Do forces for ghosts.
      potential_factor = 0.5;
      dolist<true, 1>(neighbor_list[2]);
      potential_factor = 1.0;
    }

  private:
    //! \brief Do the forces for a single list. 
    template<bool wrapping, int first_particle_type> inline void dolist(const vector<vector<int> >& neighborlist) const {
      // If the list is empty, return.
      if (neighborlist.empty()) return;

      // Get the data pointers.
      auto x0 = simData->X<first_particle_type>();
      auto x1 = simData->X();
      auto rd0 = simData->Sg<first_particle_type>();
      auto rd1 = simData->Sg();

      // Set up internal arrays.
      static_cast<const ForceType*>(this)->template set_types<first_particle_type, 0>();

      // Needed constants
      RealType R0, R1, rsqr, X[dims];

      // Get the bounds and boundary conditions...
      BCFlag bcs[dims];
      RealType widths[dims];
      // ... but only if wrapping.
      if (wrapping) {
        // Fill the vectors.
        for (int d=0; d<dims; ++d) {
          bcs[d] = gflow->getBC(d);
          widths[d] = gflow->getBounds().wd(d);
        }
      }

      for (int id0=0; id0<neighborlist.size(); ++id0) {
        auto &n_list = neighborlist[id0];
        if (n_list.empty()) continue;
        // Get radius of first particle.
        real *X0 = x0(id0);
        R0 = rd0(id0);
        // Go through the neighbors of id1.
        for (int j=0; j<n_list.size(); ++j) {
          int id1 = n_list[j];
          // Calculate displacement.
          subtract_vec<dims>(X0, x1(id1), X);
          // Harmonic corrections to distance.
          if (wrapping) harmonic_correction<dims>(bcs, X, widths);
          // Calculate squared distance
          rsqr = dot_vec<dims>(X, X);
          // Get radii
          R1 = rd1(id1);
          // If close, interact.
          if (rsqr < sqr((R0 + R1)*cutoff))
            static_cast<const ForceType*>(this)->force(id0, id1, R0, R1, rsqr, X);
          // Unset history.
          else if (using_history) {
            neighbor_history->in_contact(id0, j) = false;
            neighbor_history->value(id0, j) = 0.f;
          }
        }
      }
    }

    //! \brief A list of neighors of each particle.
    vector<vector<int> > neighbor_list[3];
    //! \brief Keep track of the number of interaction pairs.
    int _size = 0;

  };

}
#endif // __GENERIC_HANDLERS_HPP__GFLOW__