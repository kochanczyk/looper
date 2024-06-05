// This file is a part of Looper, an implementation of trajectory looping,
// distributed under the GNU General Public License, version 3.0 (see file 
// License.txt contained in the source code package or the copy online at
// http://www.gnu.org/licenses/gpl-3.0.txt).

#include "Settings.hpp"

// instantiate
size_t constexpr RunnerSettings::kSimulationSeedsBase;
size_t constexpr RunnerSettings::kNumOfReplicas;

// instantiate
double constexpr PresimulationSettings::kSearchRadius;
bool   const     PresimulationSettings::bUnloadUnnecessaryFrames;

// instantiate
double constexpr SimulationSettings::kShuffleTrajFrxnFront;
double constexpr SimulationSettings::kShuffleTrajFrxnBack;
size_t constexpr SimulationSettings::kShuffleTrajMinFramesFront;
size_t constexpr SimulationSettings::kShuffleTrajMinFramesBack;
size_t constexpr SimulationSettings::kNumOfFramePairsToMatch;
size_t constexpr SimulationSettings::kNumOfMatchingTrialDisplacements;
size_t constexpr SimulationSettings::kWriteToChemStateFileEverySteps;
#if CTHULHU_DUMPS_PDB
size_t constexpr SimulationSettings::kWriteToPDBFileEverySteps;
#endif
