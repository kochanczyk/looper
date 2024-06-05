// This file is a part of Looper, an implementation of trajectory looping,
// distributed under the GNU General Public License, version 3.0 (see file 
// License.txt contained in the source code package or the copy online at
// http://www.gnu.org/licenses/gpl-3.0.txt).

#ifndef LOOPER_SETTINGS_HPP
#define LOOPER_SETTINGS_HPP

#include "common.hpp"
#include "Chemistry.hpp"
#include <random>
#include <cassert>
#include <cstring>
#include <thread>

using chemistry_t = ChemistrySystem1;


struct RunnerSettings {
  size_t static constexpr kSimulationSeedsBase = 0u;
  size_t static constexpr kNumOfReplicas = 4u;
  size_t const kNumOfThreads = std::min(kNumOfReplicas,
                                        static_cast<size_t>(
                                             std::thread::hardware_concurrency()));
#if (LOOPER_DUMPS_PDB || LOOPER_CALCULATES_MSD)
  static_assert(kNumOfReplicas == 1, 
                "Only a single thread (single replica) can be running when "
                    "reporting MSD or dumping a PDB file.");
#endif
};


struct PresimulationSettings {
  double static constexpr kSearchRadius = 1.0 + 0.1;  // center-to-center distance
  bool   static const     bUnloadUnnecessaryFrames = true;
#if (LOOPER_DUMPS_PDB || LOOPER_CALCULATES_MSD)
  static_assert(bUnloadUnnecessaryFrames == false, 
                "Cannot report MSD or dump PDB file without coordinates.");
#endif
};


struct SimulationSettings {

  SimulationSettings(int s, const char* ofnb)
  : seed{s}, outputFileNameBase{ofnb}
  {
    assert(outputFileNameBase and (strlen(outputFileNameBase) > 0));
    assert(kTimeEnd > 0.);
  }

  // data
  int const seed = 0;                       // do not modify here
  const char* outputFileNameBase = nullptr; // do not modify here

  // simulation duration
  double const kTimeEnd = 1e5;

  // chemical rates and time scaling
  double const kLambda = 1.0e-2;

  // parameters for trajectory frames shuffling and matching
  //
  double static constexpr kShuffleTrajFrxnFront = 0.0;                  // <- If Frxn==0 & corresp.
  double static constexpr kShuffleTrajFrxnBack = kShuffleTrajFrxnFront; //    MinFrames==1, then a
  size_t static constexpr kShuffleTrajMinFramesFront = 1u;              //    single frame is taken.
  size_t static constexpr kShuffleTrajMinFramesBack = kShuffleTrajMinFramesFront;
  //
  size_t static constexpr kNumOfFramePairsToMatch  = 1u;  // per a single traj-to-traj join
  size_t static constexpr kNumOfMatchingTrialDisplacements = 1u;  // for a pair of frames
  //
  static_assert(kShuffleTrajMinFramesFront >= kNumOfFramePairsToMatch, "It's not advised "
                "to pick more frames for matching than available in the considered range.");

  // logging frequency
  size_t static constexpr kWriteToChemStateFileEverySteps = 100u;
#if LOOPER_DUMPS_PDB
  size_t static constexpr kWriteToPDBFileEverySteps = 100u;
  static_assert(kWriteToPDBFileEverySteps > 0, "PDB dumping steps interval must be > 0.");
#endif
};

#endif
