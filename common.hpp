// This file is a part of Looper, an implementation of trajectory looping,
// distributed under the GNU General Public License, version 3.0 (see file 
// License.txt contained in the source code package or the copy online at
// http://www.gnu.org/licenses/gpl-3.0.txt).
//
#ifndef LOOPER_COMMON_HPP
#define LOOPER_COMMON_HPP

#include <cstdint>
#include <random>
#if LOOPER_USES_BOOST
#  include <boost/random.hpp>
#endif

#define LOOPER_TERMINATES_GRACEFULLY  1

#define LOOPER_CALCULATES_MSD    0
#define LOOPER_DUMPS_PDB         0
#define LOOPER_DUMPS_PDB_SHIFTED 0

namespace looper {

using mol_index_t = uint16_t;

#if LOOPER_USES_BOOST
using random_number_generator_t = boost::random::mt19937;
#else
using random_number_generator_t = std::mt19937;
#endif

} // namespace looper

// definitions consistency check
#if LOOPER_DUMPS_PDB==0
# if LOOPER_DUMPS_PDB_SHIFTED==1
#  error "Conflicting PDB dumping options"
# endif
#endif

#endif

