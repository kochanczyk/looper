// This file is a part of Looper, an implementation of trajectory looping,
// distributed under the GNU General Public License, version 3.0 (see file 
// License.txt contained in the source code package or the copy online at
// http://www.gnu.org/licenses/gpl-3.0.txt).

#ifndef LOOPER_TRAJECTORY_HPP
#define LOOPER_TRAJECTORY_HPP

#include "common.hpp"
#include "TrajectoryStoragePolicies.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#if __GNUG__
#  include <malloc.h>
#endif
#include <dlib/optimization/max_cost_assignment.h>


inline
double dist_sq_in_pbc_box(double x1, double y1, double z1,
                          double x2, double y2, double z2,
                          double box_sz_x,  double box_sz_y,  double box_sz_z,
                          double half_box_sz_x, double half_box_sz_y, double half_box_sz_z)
{
  if      (x1 - x2 > half_box_sz_x) { x2 += box_sz_x; }
  else if (x2 - x1 > half_box_sz_x) { x1 += box_sz_x; }
  if      (y1 - y2 > half_box_sz_y) { y2 += box_sz_y; }
  else if (y2 - y1 > half_box_sz_y) { y1 += box_sz_y; }
  if      (z1 - z2 > half_box_sz_z) { z2 += box_sz_z; }
  else if (z2 - z1 > half_box_sz_z) { z1 += box_sz_z; }
  double const dx = x2 - x1,  dy = y2 - y1,  dz = z2 - z1;
  return dx*dx + dy*dy + dz*dz;
}


/**
 *  @brief This trajectory stores coordinates of molecules.
 */
template <class StoragePolicyT = TranslationAndRotation>
struct RawTrajectory {

  using storage_policy_t = StoragePolicyT;
  using coord_t = typename storage_policy_t::Coordinate;
  using frame_t = std::vector<coord_t>;

  static_assert(std::is_base_of<TrajectoryQuaternionStoragePolicy,storage_policy_t>::value,
                "Inadequate storage policy for instantiating a trajectory.");

  RawTrajectory(const char* filename, std::pair<int,int> range = {0,-1})
  {
    assert(std::string(filename).size() > 0);
    assert(range.first >= 0);
    assert((range.second == -1) or (range.second >= 0 and range.first <= range.second));

//  cout << "Info: Rotations handling: " << storage_policy_t::has_quaternion? "yes" : "no" << endl;

    std::ifstream fi(filename, std::fstream::binary);
    if (not fi) {
      std::stringstream what;  what <<"Cannot read file "<< filename <<".";
      throw std::runtime_error(what.str().c_str());
    }

    fi.read(reinterpret_cast<char*>(&n_mols), sizeof(n_mols));
    fi.read(reinterpret_cast<char*>(&delta_t), sizeof(delta_t));
    fi.read(reinterpret_cast<char*>(&box_sz_x), sizeof(box_sz_x));
    if (box_sz_x < 0) {
      fi.read(reinterpret_cast<char*>(&box_sz_y), sizeof(box_sz_y));
      fi.read(reinterpret_cast<char*>(&box_sz_z), sizeof(box_sz_z));
      box_sz_x = -box_sz_x;
      box_sz_y = -box_sz_y;
      box_sz_z = -box_sz_z;
    } else {
      box_sz_y =  box_sz_x;
      box_sz_z =  box_sz_x;
    }

    auto const after_hdr_pos = fi.tellg();
    fi.seekg(0, std::ios::end);
    auto const n_body_bytes = fi.tellg() - after_hdr_pos;
    fi.seekg(after_hdr_pos, std::ios::beg);
    auto const frame_data_sz = n_mols*7*sizeof(double);
    auto const frame_sz = sizeof(int32_t) + sizeof(double) + frame_data_sz;
    auto const n_estim_complete_frames = n_body_bytes / frame_sz;
    if ((n_body_bytes % frame_sz ) > 0) {
      fprintf(stderr,"Warning: Last frame appears incomplete. "
                     "Last complete-frame index is %lu.\n",n_estim_complete_frames-1);
    }
//  cout << "Estimated complete frames: " << n_estim_complete_frames << endl;

    if (range.second == -1) { range.second = n_estim_complete_frames-1; }

    for (auto const& bound: std::vector<std::pair<int,const char*>>
                            {{range.first,"lower"}, {range.second,"upper"}}) {
      if ((bound.first >= 0) and not (bound.first  < static_cast<int>(n_estim_complete_frames))) {
        std::stringstream what;
        what << "Trajectory range: "<< bound.second <<" bound index must be lower than the "
                "trajectory length (estimated to be "<< n_estim_complete_frames <<").";
        throw std::runtime_error(what.str().c_str());
      }
    } // for range bounds

    auto n_frames = n_estim_complete_frames;
    if (range.first >=0) {
      n_frames -= range.first;
      fi.seekg(static_cast<size_t>(fi.tellg()) + range.first*frame_sz, std::ios::beg);
      assert(fi.tellg() < after_hdr_pos + n_body_bytes);
    }
    if (range.second >= 0) {
      n_frames -= n_estim_complete_frames - range.second - 1;
    }
    assert(n_frames > 0);
    assert(n_frames <= n_estim_complete_frames);
//  cout << "Frames after clipping: " << n_frames << endl;

    try {
      frames = decltype(frames)(n_frames, frame_t(n_mols));
    } catch (std::bad_alloc& ba) {
      fprintf(stderr, "Error in %s: %s.\n", __func__, ba.what());
      throw;
    }

    for (auto frame_i=0u; frame_i<n_frames; ++frame_i) {

      int32_t step_i; fi.read(reinterpret_cast<char*>(&step_i), sizeof(step_i));
      double  cur_t;  fi.read(reinterpret_cast<char*>(&cur_t),  sizeof(cur_t));
//    cout << frame_i << " (" << step_i << ") @ " << cur_t << '\n';

      assert(   fi.tellg() + static_cast<decltype(n_body_bytes)>(n_mols*(7*sizeof(double)))
             <= after_hdr_pos + n_body_bytes );
      if (storage_policy_t::has_quaternion) {
        fi.read(reinterpret_cast<char*>(frames[frame_i].data()), n_mols*(7*sizeof(double)));
      } else {
        for (auto mol_i=0; mol_i<n_mols; ++mol_i) {
          fi.read(reinterpret_cast<char*>(&frames[frame_i][mol_i]), sizeof(coord_t));
          fi.seekg(static_cast<size_t>(fi.tellg())
                   + (7*sizeof(double) - sizeof(coord_t)));
        } //for each molecule
      } // if quat

    } // for each time frame

  } // RawTrajectory()

  RawTrajectory(RawTrajectory const&) = delete;
  RawTrajectory& operator=(RawTrajectory const&) = delete;

  inline size_t length() const { return frames.size(); }

  void deallocate(int frame_since, int frame_until)
  {
    for (int frame_i=frame_since; frame_i<=frame_until; ++frame_i) {
      frames[frame_i].clear();
      frames[frame_i].shrink_to_fit();
    }
#if __GNUG__
    malloc_trim(0);  // force returning memory to OS
#endif
  }

  // -- - - -  -  -  -   -

  int32_t n_mols;
  double  delta_t;
  double  box_sz_x, box_sz_y, box_sz_z;

  std::vector<frame_t> frames;

  // -- - - -  -  -  -   -

public:
  using assignment_annotated_t = std::tuple<std::vector<long>,double, coord_t>;

  template <class RandRealDistT, class RngT>
  assignment_annotated_t
  find_best_assignment(frame_t const& frame_ref, frame_t const& frame_probe, size_t n_trials,
      RandRealDistT& rand_shift_x, RandRealDistT& rand_shift_y, RandRealDistT& rand_shift_z, RngT& rng)
  const
  {
    assert(frame_ref.size() == frame_probe.size());
    std::vector<assignment_annotated_t> assignments;  assignments.reserve(n_trials);

    for (size_t trial_i = 0; trial_i < n_trials; ++trial_i) {

      // shift the probe frame
      const coord_t shift_xyz = { rand_shift_x(rng), rand_shift_y(rng), rand_shift_z(rng) };
      frame_t frame_probe_shifted = frame_probe;
      for (auto& coord: frame_probe_shifted) { coord += shift_xyz; }

      // wrap the reference frame according to the shift; old indices are preserved
      frame_t frame_ref_wrapped = frame_ref;
      for (auto& coord: frame_ref_wrapped) {
        if (not ((shift_xyz.x <= coord.x) and (coord.x < shift_xyz.x + box_sz_x))) { // needs a wrap
          coord.x += (shift_xyz.x > 0) ? +box_sz_x : -box_sz_x;
        }
        if (not ((shift_xyz.y <= coord.y) and (coord.y < shift_xyz.y + box_sz_y))) {  // needs a wrap
          coord.y += (shift_xyz.y > 0) ? +box_sz_y : -box_sz_y;
        }
        if (not ((shift_xyz.z <= coord.z) and (coord.z < shift_xyz.z + box_sz_z))) {  // needs a wrap
          coord.z += (shift_xyz.z > 0) ? +box_sz_z : -box_sz_z;
        }
      } // for each coord in the (wrapped) reference frame

      auto const a = find_assignment_with_rmsd_(frame_ref_wrapped, frame_probe_shifted);
      assignments.emplace_back(std::get<0>(a), std::get<1>(a), shift_xyz);

    } // for each trial

    // choose that of lowest RMSD
    std::sort(assignments.begin(), assignments.end(),
              [](auto const& a1, auto const& a2){ return std::get<1>(a1) < std::get<1>(a2); });
    assert(std::get<1>(assignments.front()) <= std::get<1>(assignments.back()));
    return assignments.front(); // of lowest RMSD
  }

private:
  using assignment_with_rmsd_t_  = std::tuple<std::vector<long>,double>;

  assignment_with_rmsd_t_ find_assignment_with_rmsd_(frame_t const& frameA, frame_t const& frameB)
  const
  {
    assert(frameA.size() == frameB.size());
    auto const n_mols = frameA.size();
    double const half_box_sz_x = 0.5*box_sz_x,
                 half_box_sz_y = 0.5*box_sz_y,
                 half_box_sz_z = 0.5*box_sz_z;
    double const max_box_sz = std::max(box_sz_x, std::max(box_sz_y, box_sz_z));
    double const max_d_sq = std::pow((1.0*max_box_sz)*std::sqrt(3),2);

    auto ds_sq = std::vector<std::vector<double>>(n_mols, std::vector<double>(n_mols));
    dlib::matrix<long> ds_sq_int(n_mols,n_mols);
    for (size_t i1=0, e1=n_mols; i1<e1; ++i1) {
      auto const& mol1 = frameA[i1];
      for (size_t i2=0, e2=n_mols; i2<e2; ++i2) {
        auto const& mol2 = frameB[i2];
        ds_sq[i1][i2] = dist_sq_in_pbc_box( mol1.x, mol1.y, mol1.z,
                                            mol2.x, mol2.y, mol2.z,
                                            box_sz_x,  box_sz_y,  box_sz_z,
                                            half_box_sz_x, half_box_sz_y, half_box_sz_z );
        auto d_sq_int = static_cast<long>(ds_sq[i1][i2]/max_d_sq*std::numeric_limits<long>::max());
        ds_sq_int(i1,i2) = -d_sq_int;
    } }

    std::vector<long> assignment = dlib::max_cost_assignment(ds_sq_int);  // Ta-dam!
    assert(assignment.size() == n_mols);

    auto s=0.;  for (int i=0, e=n_mols; i<e; ++i) { s += ds_sq[i][assignment[i]]; }
    auto const rmsd = s/n_mols;

    return assignment_with_rmsd_t_{assignment,rmsd};
  }

}; // class RawTrajectory


/**
 *  @brief This trajectory stores pairs of molecules that are in contact.
 *         In the article by Zuk et al., it is called 'the base contacts sequence'.
 */
template <class RawTrajectoryT>
struct ProxTrajectory {

  using raw_trajectory_t = RawTrajectoryT;
  using storage_policy_t = typename raw_trajectory_t::storage_policy_t;
  using raw_frame_t = typename raw_trajectory_t::frame_t;
  using proximity_t = typename storage_policy_t::Proximity;
  using proximities_t = std::vector<proximity_t>;
  using molecular_entries_t = std::vector<proximities_t>;  // size() == n_mols

  static_assert(std::is_base_of<TrajectoryQuaternionStoragePolicy,storage_policy_t>::value,
                "Inadequate storage policy for instantiating the proximity trajectory class:\n"
                "Sorry, cannot use foreign storage policies designed for other classes.");

  static_assert(storage_policy_t::has_quaternion? sizeof(proximity_t)==sizeof(neighbor_index_t):true,
                "When using rotation, proximity info should contain only a neighbor index.\n");

  ProxTrajectory(raw_trajectory_t const& trajectory, double cutoff)
  : n_mols  {trajectory.n_mols},
    delta_t {trajectory.delta_t},
    box_sz_x{trajectory.box_sz_x},
    box_sz_y{trajectory.box_sz_y},
    box_sz_z{trajectory.box_sz_z},
    length  {trajectory.length()},
    cutoff_sq{cutoff*cutoff}
  {
    assert(n_mols > 1);
    assert(n_mols <= std::numeric_limits<neighbor_index_t>::max());

    // -- preallocate
    try {
      frames = decltype(frames)(length,molecular_entries_t(n_mols,proximities_t{}));
    } catch (std::bad_alloc& ba) {
      fprintf(stderr,"Problem with memory around %s():%d %s.\n", __func__,__LINE__,ba.what());
      throw;
    }

    // --  fill
    double const half_box_sz_x = 0.5*box_sz_x,
                 half_box_sz_y = 0.5*box_sz_y,
                 half_box_sz_z = 0.5*box_sz_z;
    for (int fi=0, fe=length; fi<fe; ++fi) {
      raw_frame_t const& raw_frame = trajectory.frames[fi];
      auto& prox_frame = frames[fi];
      for (neighbor_index_t mol_j=0u; mol_j<n_mols; ++mol_j) {
        auto& prox_frame_mol_proximities = prox_frame[mol_j];
        auto const& jth_mol = raw_frame[mol_j];
        for (neighbor_index_t mol_k=mol_j+1; mol_k<n_mols; ++mol_k) {  // CONVENTION c2a4071a: do *not* calculated distances twice

          auto const& kth_mol = raw_frame[mol_k];
          auto const d_sq = dist_sq_in_pbc_box(jth_mol.x,jth_mol.y,jth_mol.z,
                                               kth_mol.x,kth_mol.y,kth_mol.z,
                                               box_sz_x,  box_sz_y,  box_sz_z,
                                               half_box_sz_x, half_box_sz_y, half_box_sz_z );
          if (d_sq < cutoff_sq) {
            prox_frame_mol_proximities.emplace_back(mol_k,d_sq);
            prox_frame[mol_k].emplace_back(mol_j,d_sq);
          } // if within range
        } // for mol_k
        prox_frame_mol_proximities.shrink_to_fit();
      } // for mol_j
      prox_frame.shrink_to_fit();
    } // for each frame

  } // ProxTrajectory()

  ProxTrajectory(ProxTrajectory const&) = delete;
  ProxTrajectory& operator=(ProxTrajectory const&) = delete;

  // -- - - -  -  -  -   -

  int32_t const& n_mols;
  double const& delta_t;
  double const& box_sz_x;
  double const& box_sz_y;
  double const& box_sz_z;

  size_t const length;
  double const cutoff_sq;

  std::vector<molecular_entries_t> frames;  // size() == trajectory.length()

}; // class ProxTrajectory

#endif
