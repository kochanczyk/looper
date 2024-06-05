// This file is a part of Looper, an implementation of trajectory looping,
// distributed under the GNU General Public License, version 3.0 (see file 
// License.txt contained in the source code package or the copy online at
// http://www.gnu.org/licenses/gpl-3.0.txt).
//
#include "common.hpp"
#include "Settings.hpp"
#include "Trajectory.hpp"
#include "ThreadPool.hpp"
#include <type_traits>
#include <csignal>
#include <cstdlib>
#include <chrono>

using raw_traj_t = RawTrajectory<traj_storage_t>;
using prox_traj_t = ProxTrajectory<raw_traj_t>;

#if LOOPER_TERMINATES_GRACEFULLY
std::mutex gSimulationMutex;
bool gSimulationStop = false;
#endif

// A thread-function
//
void simulation(std::vector<std::unique_ptr<const raw_traj_t>> const&  raw_trajs,
                std::vector<std::unique_ptr<const prox_traj_t>> const& prox_trajs,
                SimulationSettings settings)
{
  assert(not prox_trajs.empty());
  assert(raw_trajs.size() == prox_trajs.size());

  // -- statistics
  struct { size_t n_joins = 0u;
           size_t n_steps = 0u; } statistics;

  // -- random selectors
  //
  looper::random_number_generator_t rng(settings.seed);
  //
  // * random trajectory chooser:
  std::uniform_int_distribution<> traj_idx_picker(0, prox_trajs.size()-1);
  //
  // * random front/back frame number pickers (a pair per traj):
  std::vector< std::uniform_int_distribution<> > frame_front_idx_pickers, frame_back_idx_pickers;
  for (auto const& prox_traj: prox_trajs) {
    auto const L = prox_traj->length;
    using len_t = std::remove_const<decltype(L)>::type;
    assert(L > 0);
    int unsigned const from_front = std::max(  SimulationSettings::kShuffleTrajMinFramesFront,
                            static_cast<len_t>(SimulationSettings::kShuffleTrajFrxnFront * L) );
    int unsigned const from_back  = std::max(  SimulationSettings::kShuffleTrajMinFramesBack,
                            static_cast<len_t>(SimulationSettings::kShuffleTrajFrxnBack * L)  );
    assert((from_front < L) and (from_back < L));
    assert(0+from_front < L-1 - from_back);
    frame_front_idx_pickers.emplace_back(0, from_front-1);
    frame_back_idx_pickers.emplace_back(L-from_back,L-1);
//  cout <<"Traj #"<< frame_front_idx_pickers.size()-1 <<"["<< prox_traj.get() <<"] of length="<< L <<" "
//       <<"["<< 0 <<":"<< from_front-1 <<"]->>-["<< L-from_back <<":"<< L-1 <<"]" << endl;
  }

  // -- trajectory and its frame range (initial)
  //
  auto traj_i      = traj_idx_picker(rng);
//auto traj_i_prev = traj_i;  // UNUSED LINKED:7f36e6
  raw_traj_t const* rtraj = raw_trajs[traj_i].get();
  auto rtraj_prev = rtraj;
  prox_traj_t const*ptraj = prox_trajs[traj_i].get();
//auto ptraj_prev = ptraj;  // UNUSED
  size_t frame_front_i = 0;
  size_t frame_back_i = rtraj->length() - 1;
#if LOOPER_CALCULATES_MSD
  size_t prev_frame_back_i = std::numeric_limits<size_t>::max();
#endif

  // -- for matching
  using match_t = std::pair<size_t,typename raw_traj_t::assignment_annotated_t>;
  auto const n_mols = rtraj->n_mols;
  std::uniform_real_distribution<> rand_shift_x(-0.5*rtraj->box_sz_x,+0.5*rtraj->box_sz_x),
                                   rand_shift_y(-0.5*rtraj->box_sz_y,+0.5*rtraj->box_sz_y),
                                   rand_shift_z(-0.5*rtraj->box_sz_z,+0.5*rtraj->box_sz_z);
  std::vector<size_t> mol_i_orig2curr;
  std::vector<size_t> mol_i_curr2orig;
  mol_i_curr2orig.reserve(n_mols);  for (int i=0; i<n_mols; ++i) { mol_i_curr2orig.push_back(i); }
  mol_i_orig2curr.reserve(n_mols);  for (int i=0; i<n_mols; ++i) { mol_i_orig2curr.push_back(i); }

  // -- chemistry
#if 0
  // -- non-uniform initial conditions
  assert(rtraj->box_sz_x >  rtraj->box_sz_y);
  assert(rtraj->box_sz_y == rtraj->box_sz_z);
  std::vector<mol_index_t> t0_active;
  double const width_fxn = 0.1;  // INITIAL_CONDITION_PARAMETER GEOMETRIC
  for (mol_index_t mi = 0; mi < n_mols; ++mi) {
    // that
    if (rtraj->box_sz_x*(0.5 - width_fxn/2.) < rtraj->frames[0][mi].x
                                           and rtraj->frames[0][mi].x < rtraj->box_sz_x*(0.5 + width_fxn/2.)) {
      t0_active.push_back(mi);
  } }
  //
  chemistry_t chemistry(n_mols, settings.seed, settings.kLambda,
                        settings.outputFileNameBase, &t0_active);
  //
#else
  chemistry_t chemistry(n_mols, settings.seed, settings.kLambda,
                        settings.outputFileNameBase);
#endif

#if LOOPER_CALCULATES_MSD
  std::vector<size_t> mol_i_curr2prev;  mol_i_curr2prev.assign(n_mols,0);
  std::vector<raw_traj_t::coord_t> initial_positions, current_positions_unjoined;
  std::ofstream rmsd_ofs("MSD.dat");
#endif

#if LOOPER_CALCULATES_MSD || (LOOPER_DUMPS_PDB && LOOPER_DUMPS_PDB_SHIFTED)
  raw_traj_t::coord_t last_assignment_shift{0.,0.,0};
#endif

#if LOOPER_DUMPS_PDB

# if LOOPER_DUMPS_PDB_SHIFTED
  raw_traj_t::coord_t cumulative_assignment_shift{0.,0.,0};
# endif

  assert(SimulationSettings::kWriteToPDBFileEverySteps > 0);
  std::ofstream pdb_ofs("SIM.pdb");

  pdb_ofs <<"MODEL 1\n";
  for (int mol_i = 0; mol_i<n_mols; ++mol_i) {
    char pdb_line[100];// = {'\0'};
    auto mol_curr_pos = rtraj->frames[0][mol_i];
  //mol_curr_pos -= cumulative_assignment_shift;
    auto const mol_chem_state = chemistry.state(mol_i);
    auto const mol_obs_index = chemistry.observable_index(mol_chem_state);
  //float const beta_factor = 10.f*static_cast<float>(mol_i)/static_cast<float>(n_mols);
    float const beta_factor = 10.f*mol_obs_index/chemistry.observables_count();
    char elem = 'C';
    if (mol_obs_index == 0)  elem = 'C';
    if (mol_obs_index == 1)  elem = 'N';
    if (mol_obs_index == 2)  elem = 'O';
    if (mol_obs_index == 3)  elem = 'H';
    sprintf(pdb_line, "ATOM  %5d %1s%1c%2s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f       %2s\n",
       mol_i, " ", elem,"  ", "ALA" , "A" , mol_i,
       mol_curr_pos.x, mol_curr_pos.y, mol_curr_pos.z,  1.0, beta_factor, "C ");
    pdb_ofs << pdb_line;
  } // for
  pdb_ofs <<"ENDMDL\n";
#endif // LOOPER_DUMPS_PDB

  // BEGIN grand simulation loop
  double const delta_t = rtraj->delta_t;
  assert(delta_t > 0);
  double t = 0.;
  bool already_had_to_join_trajs = false;
  bool end_time_reached = false;
  while (t <= settings.kTimeEnd) {

    if (already_had_to_join_trajs) {

      // decide which trajectory will be next:
      traj_i = traj_idx_picker(rng);
      rtraj = raw_trajs[traj_i].get();
      ptraj = prox_trajs[traj_i].get();

//    cout <<"Info: [t="<< t <<"] Trajectory switch, current:"<< traj_i << endl;

      // pick several frames at the front of the next trajectory:
      std::set<size_t> next_traj_frame_front_indices;
      auto& next_traj_front_frame_picker = frame_front_idx_pickers[traj_i];
      while (next_traj_frame_front_indices.size() < SimulationSettings::kNumOfFramePairsToMatch) {
        next_traj_frame_front_indices.insert( next_traj_front_frame_picker(rng) );
      }

      // find the best match of the previous trajectory last frame and one of just-picked frames:
      std::vector<match_t> matches;
      auto const& last_frame = rtraj_prev->frames[frame_back_i];
      for (size_t candidate_front_i: next_traj_frame_front_indices) {
        auto const& candidate_next_frame = rtraj->frames[candidate_front_i];
        auto const a = rtraj->find_best_assignment(last_frame, candidate_next_frame,
                                                   SimulationSettings::kNumOfMatchingTrialDisplacements,
                                                   rand_shift_x, rand_shift_y, rand_shift_z, rng);
        matches.emplace_back(candidate_front_i, std::move(a));
      }
      std::sort(matches.begin(), matches.end(), [](match_t const& m1, match_t const& m2)
                                       { return std::get<1>(m1.second) < std::get<1>(m2.second); });
      auto const& best_match = matches.front();
//    cout << "Info: Best assignment corrected RMSD "
//         << get<1>(matches.back().second)  <<" --> "<< get<1>(matches.front().second)
//         <<" (Δ = ~"<< static_cast<int>((get<1>(matches.back().second)-get<1>(best_match.second))
//                                         /get<1>(matches.back().second) * 100) <<"%)"<< endl;

      // set frame range for the next trajectory
#if LOOPER_CALCULATES_MSD
      prev_frame_back_i = frame_back_i;
#endif
      frame_front_i = best_match.first;
      frame_back_i = frame_back_idx_pickers[traj_i](rng);

      // handle assignment
      auto const& assignment = std::get<0>(best_match.second);
      assert(assignment.size() == (size_t)n_mols);
      std::vector<size_t> mol_i_curr2orig_next(n_mols, std::numeric_limits<size_t>::max()); // marker
      std::vector<size_t> mol_i_orig2curr_next(n_mols, std::numeric_limits<size_t>::max()); // marker
      for (int assign_i=0; assign_i < n_mols; ++assign_i) {
        auto const mol_i_prev = assign_i;                    // in fig: 5
        auto const mol_i_curr = assignment[assign_i];        // in fig: 7
        auto const mol_i_orig = mol_i_curr2orig[mol_i_prev]; // in fig: 2
        mol_i_curr2orig_next[mol_i_curr] = mol_i_orig;       // in fig: @[7]=2
        mol_i_orig2curr_next[mol_i_orig] = mol_i_curr;       // in fig: @[2]=7
#if LOOPER_CALCULATES_MSD
        mol_i_curr2prev[mol_i_curr] = mol_i_prev;
#endif
      } // for assign_i

#if LOOPER_CALCULATES_MSD
      last_assignment_shift = std::get<2>(best_match.second);
#endif
#if LOOPER_DUMPS_PDB
# if LOOPER_DUMPS_PDB_SHIFTED
      last_assignment_shift = std::get<2>(best_match.second);
      cumulative_assignment_shift += last_assignment_shift;
# endif
#endif

      mol_i_curr2orig.swap(mol_i_curr2orig_next);  // efficient op=
      mol_i_orig2curr.swap(mol_i_orig2curr_next);  // efficient op=
      assert(!std::count(mol_i_curr2orig.begin(),mol_i_curr2orig.end(),
             std::numeric_limits<size_t>::max()));
      assert(!std::count(mol_i_orig2curr.begin(),mol_i_orig2curr.end(),
              std::numeric_limits<size_t>::max()));

      statistics.n_joins++;
    } // if already joined trajectories

    // BEGIN continuous trajectory processing
    for (auto frame_i = frame_front_i;  frame_i <= frame_back_i;  ++frame_i) {
      auto const& frame = ptraj->frames[frame_i];

#if LOOPER_CALCULATES_MSD
      std::vector<double> cur_frame_sds;
#endif
#if LOOPER_DUMPS_PDB
      pdb_ofs <<"MODEL "<< 2+statistics.n_steps/SimulationSettings::kWriteToPDBFileEverySteps <<"\n";
#endif

      chemistry.initialize_timepoint(t);

      for (int mol_i=0; mol_i<n_mols; ++mol_i) { // ---
        auto const& mol_proximities = frame[mol_i];
//      auto const& mol_coords = rtraj->frames[frame_i][mol_i];
//      printf("%8.6f  [traj %d]  [frame %lu]  [molecule %3d]  (%6.3f %6.3f %6.3f)  {neighs:%2lu}\n",
//              t,traj_i,frame_i,mol_i,mol_coords.x,mol_coords.y,mol_coords.z,mol_proximities.size());

        // -- chemistry
        for (auto proxi: mol_proximities) {
          if (mol_i > proxi.idx) { continue; } // c2a4071a BACK
#ifndef NDEBUG
          if (not PresimulationSettings::bUnloadUnnecessaryFrames) {
            auto const& mol_coords   = rtraj->frames[frame_i][mol_i];
            auto const& neigh_coords = rtraj->frames[frame_i][proxi.idx];
            double const box_sz_x = rtraj->box_sz_x,
                         box_sz_y = rtraj->box_sz_y,
                         box_sz_z = rtraj->box_sz_z;
            auto const dsq = dist_sq_in_pbc_box(mol_coords.x, mol_coords.y, mol_coords.z,
                                                neigh_coords.x,neigh_coords.y,neigh_coords.z,
                                                box_sz_x,    box_sz_y,    box_sz_z,
                                                box_sz_x/2., box_sz_y/2., box_sz_z/2.);
            assert(proxi.dist == std::sqrt(dsq));
          }
#endif // !NDEBUG
          auto const mi1 = mol_i_curr2orig[mol_i];
          auto const mi2 = mol_i_curr2orig[proxi.idx];
          chemistry(mi1, mi2);
        } // for each molecule in proximity

#if LOOPER_CALCULATES_MSD
        if (initial_positions.size() < static_cast<size_t>(n_mols)) {
          assert((statistics.n_joins == 0) and (frame_i == 0));
          auto const& mol_curr_pos = rtraj->frames[frame_i][mol_i];
          initial_positions.push_back(raw_traj_t::coord_t{mol_curr_pos.x, mol_curr_pos.y, mol_curr_pos.z});
          current_positions_unjoined.push_back(initial_positions.back());
        } else {
          assert(initial_positions.size() == static_cast<size_t>(n_mols));

          raw_traj_t::coord_t mol_curr_pos = rtraj->frames[frame_i][mol_i];
          raw_traj_t::coord_t mol_prev_pos;
          if (frame_i == 0) {
            assert(statistics.n_joins > 0);
            mol_curr_pos += last_assignment_shift;
            mol_prev_pos =  rtraj_prev->frames[prev_frame_back_i][ mol_i_curr2prev[mol_i] ];
          } else {
            assert(frame_i > 0);
            mol_prev_pos = rtraj->frames[frame_i-1][mol_i];
          } // if

          raw_traj_t::coord_t d = mol_curr_pos;  d -= mol_prev_pos;
          if      (d.x > +rtraj->box_sz_x/2.) { d.x -= rtraj->box_sz_x; }
          else if (d.x < -rtraj->box_sz_x/2.) { d.x += rtraj->box_sz_x; }
          if      (d.y > +rtraj->box_sz_y/2.) { d.y -= rtraj->box_sz_y; }
          else if (d.y < -rtraj->box_sz_y/2.) { d.y += rtraj->box_sz_y; }
          if      (d.z > +rtraj->box_sz_z/2.) { d.z -= rtraj->box_sz_z; }
          else if (d.z < -rtraj->box_sz_z/2.) { d.z += rtraj->box_sz_z; }

          auto mol_i_orig = mol_i_curr2orig[mol_i];
          auto& mol_unwr_pos = current_positions_unjoined[mol_i_orig];
          mol_unwr_pos += d;
          auto const& mol_init_pos = initial_positions[mol_i_orig];
          raw_traj_t::coord_t D = mol_unwr_pos;  D -= mol_init_pos;
          auto const sd = D.x*D.x + D.y*D.y + D.z*D.z;
          cur_frame_sds.push_back(sd);
        } // if
#endif // LOOPER_CALCULATES_MSD

#if LOOPER_DUMPS_PDB
        if (not (statistics.n_steps % SimulationSettings::kWriteToPDBFileEverySteps)) {
          char pdb_line[100];// = {'\0'};
          auto mol_curr_pos = rtraj->frames[frame_i][ mol_i_curr2orig[mol_i] ];
# if LOOPER_DUMPS_PDB_SHIFTED
          mol_curr_pos -= cumulative_assignment_shift;
          auto const mol_chem_state = chemistry.state(                 mol_i  );
# else
          auto const mol_chem_state = chemistry.state( mol_i_curr2orig[mol_i] );
# endif
          auto const mol_obs_index = chemistry.observable_index(mol_chem_state);
        //float const beta_factor = 10.f*static_cast<float>(mol_i)/static_cast<float>(n_mols);
          float const beta_factor = 10.f*mol_obs_index/chemistry.observables_count();
          char elem = 'C';
          if (mol_obs_index == 0)  elem = 'C';
          if (mol_obs_index == 1)  elem = 'N';
          if (mol_obs_index == 2)  elem = 'O';
          if (mol_obs_index == 3)  elem = 'H';
          sprintf(pdb_line, "ATOM  %5d %1s%1c%2s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f       %2s\n",
             mol_i, " ", elem,"  ", "ALA" , "A" , mol_i,
             mol_curr_pos.x, mol_curr_pos.y, mol_curr_pos.z,  1.0, beta_factor, "C ");
          pdb_ofs << pdb_line;
        }
#endif // LOOPER_DUMPS_PDB

      } // for each molecule ---

      //prox_traj_t::molecular_entries_t const& P = frame;
      //std::vector<std::vector<traj_storage_t::Proximity>> const& PP = frame;

      chemistry.finalize_timepoint(delta_t, &frame);

#if LOOPER_CALCULATES_MSD
      assert(statistics.n_steps ? cur_frame_sds.size() == static_cast<size_t>(n_mols) : true);
      double sum_frame_sds = 0.;
      for (auto const& sd: cur_frame_sds) { sum_frame_sds += sd; }
      double const frame_msd = sum_frame_sds / n_mols;
      rmsd_ofs << t << '\t' << (frame_msd) << std::endl;
#endif // LOOPER_CALCULATES_MSD

#if LOOPER_DUMPS_PDB
      if (not (statistics.n_steps % SimulationSettings::kWriteToPDBFileEverySteps)) {
        pdb_ofs << "ENDMDL\n";
      }
#endif // LOOPER_DUMPS_PDB

      // -- save chemical state (before advancing time)
      if (not (statistics.n_steps % SimulationSettings::kWriteToChemStateFileEverySteps)) {
        chemistry.append_to_file(t);
      }

      if (t > settings.kTimeEnd) {
        end_time_reached = true;
        break;
      } else {
        t += delta_t;
        statistics.n_steps++;
      }

#if LOOPER_TERMINATES_GRACEFULLY      
      {
        std::unique_lock<std::mutex> __{gSimulationMutex};
        if (gSimulationStop) { end_time_reached = true; }
      }
#endif

    } // for each frame
    // END continuous trajectory processing

    if (end_time_reached) { break; }

  //traj_i_prev = traj_i; // UNUSED LINKED:7f36e6
    rtraj_prev  = rtraj;
  //ptraj_prev  = ptraj;  // UNUSED
    already_had_to_join_trajs  = true;
  }
  // END grand simulation loop

#if LOOPER_DUMPS_PDB
  pdb_ofs << "END\n";
#endif // LOOPER_DUMPS_PDB

  // statistics write-out:
  std::stringstream stat_ofn;
  stat_ofn << settings.outputFileNameBase <<"-sim_"<< settings.seed <<"-stats.log";
  std::ofstream(stat_ofn.str().c_str()) << "n_steps: " << statistics.n_steps << std::endl
                                        << "n_joins: " << statistics.n_joins << std::endl;
}


#if LOOPER_TERMINATES_GRACEFULLY
static char goodbye_status[256];

void use_terminator(int)
{
  static bool already_called = false;
  printf("\n[TERMINATING...%s]\n", 
         already_called ? " (PATIENCE, I WILL STOP AT THE NEXT JOIN)" : "");  
  fflush(stdout);
  {
    std::unique_lock<std::mutex> __{gSimulationMutex};
    gSimulationStop = true;
  }
  already_called = true;
}
#endif


void show_usage(char const* argv0)
{
  printf("\nUsage:\n");
  printf("      %s trajectory_file chem_prefix\n\n", argv0);
}

int main(int argc, char* argv[])
{
  // -- Parse command line arguments to obtain i/o file names:
  //
  if (argc <= 2) {
    show_usage(argv[0]);
    exit(2);
  }
  //
  std::vector<const char*> traj_input_files_names;
  for (int ai=1; ai<argc - 1; ++ai) { traj_input_files_names.push_back(argv[ai]); }
  //
  const char* output_file_name_base = argv[argc - 1];

  // -- Load raw trajectories and derive proximity trajectories:
  //
  std::vector<std::unique_ptr<const raw_traj_t>> raw_trajs;
  std::vector<std::unique_ptr<const prox_traj_t>> prox_trajs;
  for (auto const& fn: traj_input_files_names) {

    // -- Load raw trajectory frames from file
    //
    printf("Info: Loading %s ...\n", fn);  fflush(stdout);
    raw_trajs.emplace_back(new raw_traj_t{fn});
    printf("Info: Loading %s ::: done\n", fn);  fflush(stdout);

    // -- Perform several consistency checks
    //
    if (raw_trajs.back()->length() <  SimulationSettings::kShuffleTrajMinFramesFront 
                                    + SimulationSettings::kShuffleTrajMinFramesBack) {
      throw std::logic_error("Trajectory too short: the number of frames available for "
                             "match-looping is below the minimal expected frame number.");
    }
    if (raw_trajs.size() >= 2) {
      auto const& curr_traj =    *raw_trajs.crbegin();
      auto const& prev_traj = *(++raw_trajs.crbegin());
      if (curr_traj->delta_t != prev_traj->delta_t) 
      { throw std::logic_error("Mismatch of delta_t between the just-loaded & the previous traj.");}
      if (curr_traj->n_mols != prev_traj->n_mols) 
      { throw std::logic_error("Mismatch of n_mols between the just-loaded & the previous traj."); }
      if (curr_traj->box_sz_x != prev_traj->box_sz_x)
      { throw std::logic_error("Mismatch of box_sz_x between the just-loaded & the previous traj.");}
      if (curr_traj->box_sz_y != prev_traj->box_sz_y)
      { throw std::logic_error("Mismatch of box_sz_y between the just-loaded & the previous traj.");}
      if (curr_traj->box_sz_z != prev_traj->box_sz_z)
      { throw std::logic_error("Mismatch of box_sz_z between the just-loaded & the previous traj.");}
    }

    // -- Derive proximities/"proxy" trajectory
    //
    printf("Info: Generating contact sequence (r = %.2f) ...\n",PresimulationSettings::kSearchRadius);
    fflush(stdout);
    prox_trajs.emplace_back(new prox_traj_t{*raw_trajs.back(),PresimulationSettings::kSearchRadius});
    printf("Info: Generating contact sequence (r = %.2f) ::: done\n",PresimulationSettings::kSearchRadius);
    fflush(stdout);

    // -- (optional) Forget raw trajectory frames not necessary for matching frames upon rewind
    //
    if (PresimulationSettings::bUnloadUnnecessaryFrames) {
      if (    SimulationSettings::kShuffleTrajFrxnFront == 0.0
          and SimulationSettings::kShuffleTrajFrxnBack  == 0.0
          and SimulationSettings::kShuffleTrajMinFramesFront == 1u
          and SimulationSettings::kShuffleTrajMinFramesBack  == 1u) {
        auto p_last_raw_traj = const_cast<raw_traj_t*>(raw_trajs.back().get()); // CONST-CORRECTNESS VIOLATION
        int const front_frames_left = SimulationSettings::kShuffleTrajMinFramesFront,
                  back_frames_left  = SimulationSettings::kShuffleTrajMinFramesBack;
        printf("Info: Unloading frames (indices 0+%d:%zu-%d) ...\n",
                front_frames_left, p_last_raw_traj->length(), back_frames_left);   fflush(stdout);
        p_last_raw_traj->deallocate(front_frames_left, p_last_raw_traj->length()-1-back_frames_left);
        printf("Info: Unloading frames (indices 0+%d:%zu-%d) ::: done\n",
                front_frames_left, p_last_raw_traj->length(), back_frames_left);   fflush(stdout);
      } // if specific end-to-front matching conditions met
    } // if UnloadUnnecessary

  } // for each trajectory file name

#if LOOPER_TERMINATES_GRACEFULLY
  // -- Set signal handlers
  signal(SIGINT,  use_terminator);
  signal(SIGTERM, use_terminator);
#endif

  // -- Spawn several simulations running in parallel
  //
  RunnerSettings run_setts;
  ThreadPool pool{run_setts.kNumOfThreads};
  std::vector<std::future<void>> future_results;

  const char* const env_seed_base = std::getenv("LOOPER_SEED_BASE");
  int sim_seed_begin =  env_seed_base ? atoi(env_seed_base)
                                      : RunnerSettings::kSimulationSeedsBase;
  for (int sim_seed     = sim_seed_begin, 
           sim_seed_end = sim_seed_begin + RunnerSettings::kNumOfReplicas; 
           sim_seed < sim_seed_end; ++sim_seed) {
    SimulationSettings sim_setts{sim_seed, output_file_name_base};
    future_results.emplace_back( 
        pool.enqueue(simulation, std::cref(raw_trajs), std::cref(prox_trajs), sim_setts) );
  }

  printf("Info: Running %zu simulation%s on %zu thread%s ...\n",
          RunnerSettings::kNumOfReplicas, RunnerSettings::kNumOfReplicas > 1 ? "s":"",
          run_setts.kNumOfThreads, run_setts.kNumOfThreads > 1 ? "s":"");  fflush(stdout);
#if LOOPER_TERMINATES_GRACEFULLY
  sprintf(goodbye_status,
          "Info: Running %zu simulation%s on %zu thread%s ::: DONE\n",
          RunnerSettings::kNumOfReplicas, RunnerSettings::kNumOfReplicas>1 ? "s":"",
          run_setts.kNumOfThreads, run_setts.kNumOfThreads > 1 ? "s":""            );
  std::atexit([](){ printf("%s",goodbye_status);  fflush(stdout); });
#endif

  for (auto& r: future_results) { r.get(); }

  return EXIT_SUCCESS;
}

