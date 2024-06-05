// This file is a part of Looper, an implementation of trajectory looping,
// distributed under the GNU General Public License, version 3.0 (see file 
// License.txt contained in the source code package or the copy online at
// http://www.gnu.org/licenses/gpl-3.0.txt).

#ifndef LOOPER_TRAJECTORY_STORAGE_POLICIES_HPP
#define LOOPER_TRAJECTORY_STORAGE_POLICIES_HPP

#include "common.hpp"

using looper::mol_index_t;
using neighbor_index_t = mol_index_t;

struct TrajectoryQuaternionStoragePolicy { };

struct TranslationOnly : TrajectoryQuaternionStoragePolicy
{
  enum : bool { has_quaternion = false };
  struct Coordinate { double x,y,z;  // has no quaternion
    Coordinate& operator+=(Coordinate const& c) {  x += c.x;  y += c.y;  z += c.z;  return *this;  }
    Coordinate& operator-=(Coordinate const& c) {  x -= c.x;  y -= c.y;  z -= c.z;  return *this;  }
  };
  struct Proximity  {
    inline Proximity(neighbor_index_t i, double d_sq) : idx(i), dist(std::sqrt(d_sq)) { }
    neighbor_index_t const idx;  double const dist;
  };
};

struct TranslationAndRotation : TrajectoryQuaternionStoragePolicy
{
  enum : bool { has_quaternion = true };
  struct Coordinate { double x,y,z, q1=0.,q2=0.,q3=0.,q4=0.; };
  struct Proximity  {
    inline Proximity(neighbor_index_t i, double) noexcept : idx(i) { }
    neighbor_index_t const idx;  // Proximity in TranslationAndRotation has no distance.
  };
};

using traj_storage_t = TranslationOnly;

#endif
