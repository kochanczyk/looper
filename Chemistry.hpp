// This file is a part of Looper, an implementation of trajectory looping,
// distributed under the GNU General Public License, version 3.0 (see file 
// License.txt contained in the source code package or the copy online at
// http://www.gnu.org/licenses/gpl-3.0.txt).
//
#ifndef LOOPER_CHEMISTRY_HPP
#define LOOPER_CHEMISTRY_HPP

#include "common.hpp"
#include "TrajectoryStoragePolicies.hpp"
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <cmath>
#include <memory>
#include <cassert>
#include <sstream>
#include <fstream>

using looper::mol_index_t;

template <class LabelT = unsigned int>
struct ChemicalState : public std::vector<LabelT> {
  using label_t = LabelT;
  using std::vector<label_t>::vector;
};

template <class ChemicalStateT = ChemicalState<>>
struct Reagent : public std::tuple<typename ChemicalStateT::label_t*, mol_index_t,
                                   typename ChemicalStateT::label_t > {
  using label_t = typename ChemicalStateT::label_t;
  using parent_t = std::tuple<label_t*, mol_index_t, label_t>;
  using typename parent_t::tuple;

#if ((__GNUC__ >= 6) || defined(__clang__))
  Reagent(label_t* mod, mol_index_t idx, label_t ref) : parent_t{mod,idx,ref} { }
  Reagent(Reagent&& other) = default;
  Reagent(Reagent const& other) = default;
  Reagent& operator=(Reagent const& other) = default;
  Reagent& operator=(Reagent&& other) = default;
#endif

  inline label_t const* state()   const { return std::get<label_t*>(*this);    }
  inline label_t*       state()         { return std::get<label_t*>(*this);    }
  inline mol_index_t    index()   const { return std::get<mol_index_t>(*this); }
  inline label_t        stencil() const { return std::get<label_t>(*this);     }  
  inline bool           intact()  const { return *state() == stencil();        }
};


template <class ChemicalStateT>
struct Chemistry {

  using chemical_state_t = ChemicalStateT;
  using label_t = typename chemical_state_t::label_t;
  using observables_t = std::vector<label_t>;
  using mol_entries_t = std::vector<std::vector<traj_storage_t::Proximity>>;

  Chemistry(Chemistry const&)            = delete;
  Chemistry& operator=(Chemistry const&) = delete;
  virtual ~Chemistry() { }

  template<class LabelT, class... LabelsT>
  void observe(LabelT obs_head, LabelsT... obs_tail)
  {  observables_.push_back(obs_head);  observe(obs_tail...);  }
  //
  void observe() {}

  inline int observable_index(label_t obs) const  // PERFORMANCE NOTE: linear COMPLEXITY
  {
    return std::distance( observables_.cbegin(),
                          std::find(observables_.cbegin(),observables_.cend(),obs) );
  }

  inline int observables_count() const { return observables_.size(); }

  inline label_t state(mol_index_t mi) const { return chem_state_[mi]; }

  void chem_state_replace( looper::random_number_generator_t& rng,
                           label_t to_replace, size_t amount, label_t state )
  {
    assert(std::count(chem_state_.cbegin(),chem_state_.cend(), to_replace) >= (long)amount);
    auto static rand_mol_index_gen = std::uniform_int_distribution<>(0, nmols_ - 1);
    while (amount > 0) {
      auto const i = rand_mol_index_gen(rng);
      if (chem_state_[i] == to_replace) {
        chem_state_[i] = state;
        --amount;
    } }
  }

  virtual void initialize_timepoint(double t) = 0;
  virtual void finalize_timepoint(double delta_t, mol_entries_t const*) = 0;

  virtual void append_to_file(double t) const
  {
    if (chem_ofs_) {
      assert(not observables_.empty());
      *chem_ofs_ << std::scientific << lambda_*t << " ";
      for (auto obs: observables_) {
        auto const cnt = std::count(chem_state_.cbegin(), chem_state_.cend(), obs);
        *chem_ofs_ << cnt << ' ';
      }
      *chem_ofs_ << std::endl;
    }
  }

protected:
  observables_t observables_;
  chemical_state_t chem_state_;
  std::unique_ptr<std::ofstream> chem_ofs_ = nullptr;
  size_t const nmols_;
  double t_ = -1.0;
  double lambda_ = 1e-2;

  Chemistry(size_t nmols, unsigned int run_id, double lambda, 
            const char* chem_ofn_base = nullptr)
  : nmols_{nmols},
    lambda_{lambda}
  {
    assert(nmols > 0);
    assert(lambda > 0.);
    chem_state_.reserve(nmols);
    if (chem_ofn_base) {
      std::stringstream chem_ofn;
      chem_ofn << chem_ofn_base <<"-sim_"<< run_id <<"-chem.dat";
      chem_ofs_.reset(new std::ofstream (chem_ofn.str().c_str()));
    }
  }

}; // Chemistry



template <class ChemicalStateT>
struct DeterministicChemistry : public Chemistry<ChemicalStateT> {

  using parent_chemistry_t = Chemistry<ChemicalStateT>;
  using chemical_state_t = ChemicalStateT;
  using label_t = typename chemical_state_t::label_t;
  using mol_entries_t = typename parent_chemistry_t::mol_entries_t;

  DeterministicChemistry(DeterministicChemistry const&) = delete;
  DeterministicChemistry& operator=(DeterministicChemistry const&) = delete;
  virtual ~DeterministicChemistry() { }

  void initialize_timepoint(double t) { parent_chemistry_t::t_ = t; }
  void finalize_timepoint(double /*delta_t*/, mol_entries_t const* = nullptr) {}

  virtual void operator()(mol_index_t) {};                //  \__ These 2 operators should
  virtual void operator()(mol_index_t,mol_index_t) = 0;   //  /   immediately fire reaction(s).

protected:
  DeterministicChemistry(size_t nmols, unsigned int run_id, double lambda,
                         const char* chem_ofn_base = nullptr)
  : parent_chemistry_t(nmols, run_id, lambda, chem_ofn_base)
  { }

}; // DeterministicChemistry



template <class ChemicalStateT, class RandomNumGenT>
struct StochasticChemistry : public Chemistry<ChemicalStateT> {

  using chemical_state_t = ChemicalStateT;
  using parent_chemistry_t = Chemistry<chemical_state_t>;
  using label_t = typename parent_chemistry_t::label_t;
  using mol_entries_t = typename parent_chemistry_t::mol_entries_t;
  using rng_t = RandomNumGenT;

  using reagent_t = Reagent<chemical_state_t>;
  using action_t = std::function<void(reagent_t&, reagent_t&)>;
  using check_t  = std::function<bool(reagent_t const&, reagent_t const&)>;
  using event_t = std::tuple<reagent_t, reagent_t, std::pair<action_t,double>/*, check_t*/>;

  StochasticChemistry(StochasticChemistry const&) = delete;
  StochasticChemistry& operator=(StochasticChemistry const&) = delete;
  virtual ~StochasticChemistry() { }

  void initialize_timepoint(double t)
  {
    parent_chemistry_t::t_ = t;
    events_.reserve(0.25 * parent_chemistry_t::nmols_);  // PERFORMANCE HEURISTIC PARAMETER 977d44
    assert(events_.empty());
  }

  /// @brief Gillespie-type algorithm
  void finalize_timepoint(double delta_t, mol_entries_t const* neighbors = nullptr)
  {
    assert(neighbors);

    double t = 0.;
    double sum_rates = std::accumulate( events_.begin(), events_.end(), 0., 
            [this](double partial_sum, event_t const& ev) { return partial_sum + rate(ev); } );

    while (t <= delta_t) { // Gillespie algorithm over delta_t

      if (events_.empty()) { break; }
 
      assert(sum_rates > 0.);
      double const dt = 1./sum_rates * log(1./uniform01_(rng_));
      if (t + dt > delta_t) { break; }  // early break

      // -- select an event
      double const rho = sum_rates * uniform01_(rng_);
      double partial_sum_rates = 0.;
      auto const event_it = std::find_if(events_.begin(), events_.end(),
          [&partial_sum_rates,rho,this](auto const& event)
          { return (partial_sum_rates += this->rate(event)) >= rho; });
      if (event_it == events_.cend()) {
        sum_rates = std::accumulate( events_.begin(), events_.end(), 0.,
                    [this](double part_sum, event_t const& ev) { return part_sum + rate(ev); } );
        continue;
      }

      // -- increase time (here, because rates might have been recalc'ed with continue above)
      t += dt;

      // -- fire an event
      assert(event_it != events_.cend());
      event_t& event = *event_it;
      assert(applicable(event));
      if (!applicable(event)) { printf("!"); fflush(stdout); }
      execute(event);

      // -- remove all events involving the changed substrate
      auto const substrate_index = reagent2(event).index();
      for (auto evi = events_.begin(); evi != events_.end(); /* */) {
        auto& event = *evi;
        if (   reagent1(event).index() == substrate_index 
            or reagent2(event).index() == substrate_index) {
          sum_rates -= rate(event);
          std::swap(event, events_.back());
          events_.pop_back();
        } else {
          ++evi;
        }
      } // for

      // -- subscribe new events that involve the changed substrate
      auto const& prox_neighbors = (*neighbors)[substrate_index];
      for (auto prox_neighbor: prox_neighbors) {
        sum_rates += match_events(substrate_index,prox_neighbor.idx);
        sum_rates += match_events(prox_neighbor.idx,substrate_index);
      }

    } // while t < delta_t
    
    events_.clear();

  } // finalize_timepoint()

  virtual void operator()(mol_index_t) { }                 //  \__ These 2 operators should create
  virtual void operator()(mol_index_t miA,mol_index_t miB) //  /   deferred events to be fired in
  {                                                        //      finalize_timepoint().
    match_events(miA, miB);
    match_events(miB, miA);
  }

  virtual double match_events(mol_index_t mi1, mol_index_t mi2) = 0;
    

protected:
  StochasticChemistry(size_t nmols, unsigned int run_id, double lambda,
                      const char* chem_ofn_base = nullptr)
  : parent_chemistry_t{nmols, run_id, lambda, chem_ofn_base},
    rng_{run_id},
    uniform01_{0.,1.}
  { }

  inline double add_event(void (*axn)(reagent_t&,reagent_t&),
                          mol_index_t mi1, mol_index_t mi2,  double rate)
  {
    using std::placeholders::_1;
    using std::placeholders::_2;
    label_t& s1 = parent_chemistry_t::chem_state_[mi1];
    label_t& s2 = parent_chemistry_t::chem_state_[mi2];
    events_.emplace_back( std::move(reagent_t{&s1, mi1, label_t{s1}}),
                          std::move(reagent_t{&s2, mi2, label_t{s2}}),
                          std::make_pair(std::bind(axn,_1,_2), rate) );
    return rate;
  }

  rng_t rng_;
  std::uniform_real_distribution<> uniform01_;
  std::vector<event_t> events_;

  inline reagent_t const& reagent1(event_t const& e) const { return std::get<0>(e);        }
  inline reagent_t&       reagent1(event_t&       e)       { return std::get<0>(e);        }
  inline reagent_t const& reagent2(event_t const& e) const { return std::get<1>(e);        }
  inline reagent_t&       reagent2(event_t&       e)       { return std::get<1>(e);        }
  inline auto               action(event_t&       e) const { return std::get<2>(e).first;  }
  inline auto                 rate(event_t const& e) const { return std::get<2>(e).second; }
  inline bool applicable(event_t const&e)const{return reagent1(e).intact() and reagent2(e).intact();}
  inline void    execute(event_t&      e)     { action(e)(reagent1(e),reagent2(e)); }

}; // StochasticChemistry


// =================================================================================================


struct NoChemistry : public DeterministicChemistry<ChemicalState<>> {

  using parent_chemistry_t = DeterministicChemistry<ChemicalState<>>;

  NoChemistry(size_t nmols, unsigned int run_id, double lambda,
              const char* chem_ofn_base = nullptr)
  : parent_chemistry_t{nmols, run_id, lambda, chem_ofn_base}
  { }

  void operator()(mol_index_t /*mi1*/, mol_index_t /*mi2*/)
  { }

}; // NoChemistry



struct ChemistryMidasTouch : public DeterministicChemistry<ChemicalState<>> {

  using parent_chemistry_t = DeterministicChemistry<ChemicalState<>>;

  ChemistryMidasTouch(size_t nmols, unsigned int run_id, double lambda,
                      const char* chem_ofn_base = nullptr)
  : parent_chemistry_t{nmols, run_id, lambda, chem_ofn_base}
  {
    observe(0,1,2);
    label_t const bulk = static_cast<label_t>(0), distinct = static_cast<label_t>(1);
    chem_state_.assign(nmols,bulk);
    looper::random_number_generator_t rng{run_id};
    chem_state_replace(rng, bulk, 1, distinct);
  }

  void operator()(mol_index_t mi1, mol_index_t mi2)
  {
    label_t& S1 = chem_state_[mi1];
    label_t& S2 = chem_state_[mi2];
    if (S1 == 1) { S2 = 2; }
    if (S2 == 1) { S1 = 2; }
  }

}; // ChemistryMidasTouch



struct ChemistryForNewEncounters : public DeterministicChemistry<ChemicalState<>> {

  using parent_chemistry_t = DeterministicChemistry<chemical_state_t>;

  ChemistryForNewEncounters(size_t nmols, unsigned int run_id,  double lambda,
                            const char* chem_ofn_base = nullptr)
  : parent_chemistry_t{nmols, run_id, lambda, chem_ofn_base}
  {
    for (size_t i = 0; i < nmols; ++i) { observables_.push_back(static_cast<label_t>(i)); }

    label_t const initial_bulk = static_cast<label_t>(0);
    chem_state_.assign(nmols,initial_bulk);

    for (size_t i = 0; i < nmols; ++i) {
      react_times_.emplace_back(std::vector<double>(i+1, -1.0));
    }
    //
    if (chem_ofn_base) {
      std::stringstream cont_ofn;
      cont_ofn << chem_ofn_base <<"-sim_"<< run_id <<"-trea.dat";
      react_times_ofs_.reset(new std::ofstream{cont_ofn.str().c_str()});
    }
  }

  ~ChemistryForNewEncounters()
  {
    if (react_times_ofs_) {
      for (size_t i=0; i<nmols_; i++) {
        for (size_t j=0; j<=i; j++) {
          if ( react_times_[i][j] >= 0. ) { *react_times_ofs_ << react_times_[i][j] << '\n'; }
      } } // for,for
    } // if
    react_times_ofs_->flush();
  }

  void operator()(mol_index_t mi1, mol_index_t mi2)
  {
    double const react_time_mi1mi2 = (mi1 > mi2) ? react_times_[mi1][mi2]
                                                 : react_times_[mi2][mi1];
    if (react_time_mi1mi2 < 0.) {
      label_t& S1 = chem_state_[mi1];
      label_t& S2 = chem_state_[mi2];
      ++S1;
      ++S2;
      if (mi1 > mi2) { react_times_[mi1][mi2] = t_; }
      else           { react_times_[mi2][mi1] = t_; }
    }
  }

protected:
  std::unique_ptr<std::ofstream> react_times_ofs_ = nullptr;
  std::vector<std::vector<double>> react_times_;

}; // ChemistryForNewEncounters



enum struct MolStateSystem1 : unsigned char { Su, Sp, K, P };
//
struct ChemistrySystem1 : public StochasticChemistry< ChemicalState<MolStateSystem1>,
                                                      looper::random_number_generator_t > {
  using parent_chemistry_t = StochasticChemistry< ChemicalState<MolStateSystem1>,
                                                  looper::random_number_generator_t >;
  using chemical_state_t = typename parent_chemistry_t::chemical_state_t;
  using L = typename chemical_state_t::label_t;

  double const rate_q0;
  double const rate_k0;

  ChemistrySystem1(size_t nmols, unsigned int run_id, double lambda,
                   const char* chem_ofn_base = nullptr)
  : parent_chemistry_t{nmols, run_id, lambda, chem_ofn_base},
    rate_q0(lambda),
    rate_k0(0.01*rate_q0)
  {
    observe(L::Su, L::Sp, L::K, L::P);
    chem_state_.assign(nmols,L::Sp);
    chem_state_replace(rng_, L::Sp, static_cast<size_t>(0.3 *nmols), L::K);
    chem_state_replace(rng_, L::Sp, static_cast<size_t>(0.03*nmols), L::P);
  }

protected:

  double match_events(mol_index_t mi1, mol_index_t mi2)
  {
    L const s1 = chem_state_[mi1],  s2 = chem_state_[mi2];
    if      (s1==L::K and s2==L::Su) { return add_event(  phospho, mi1, mi2, rate_k0); }
    else if (s1==L::P and s2==L::Sp) { return add_event(dephospho, mi1, mi2, rate_q0); }
    else { return 0.; }
  }

  static void   phospho(reagent_t& K, reagent_t& S)
  {  assert(*K.state()==L::K and *S.state()==L::Su);  *S.state()=L::Sp;  (void)(K);  }

  static void dephospho(reagent_t& P, reagent_t& S)
  {  assert(*P.state()==L::P and *S.state()==L::Sp);  *S.state()=L::Su;  (void)(P);  }

}; // ChemistrySystem1



enum struct MolStateSystem2 : unsigned char { Ku, Kp, Kpp, P, C };
//
struct ChemistrySystem2 : public StochasticChemistry< ChemicalState<MolStateSystem2>,
                                                      looper::random_number_generator_t > {
  using parent_chemistry_t =StochasticChemistry< ChemicalState<MolStateSystem2>,
                                                 looper::random_number_generator_t >;
  using chemical_state_t = typename parent_chemistry_t::chemical_state_t;
  using L = typename chemical_state_t::label_t;

  double const rhoP, rhoK;
  double const rate_q0, rate_c1, rate_c2, rate_c3;

  ChemistrySystem2(size_t nmols, unsigned int run_id, double lambda,
                   const char* chem_ofn_base = nullptr,
                   std::vector<mol_index_t>* t0_active = nullptr) // special initial conditions
  : parent_chemistry_t{nmols, run_id, lambda, chem_ofn_base},
    rhoP{0.2}, rhoK{1.0 - rhoP},
    rate_q0{lambda}, rate_c1{0.025*rate_q0}, rate_c2{0.0075*rate_q0}, rate_c3{rate_q0}
  {
    observe(L::Ku, L::Kp, L::Kpp, L::P);

    chem_state_.assign(nmols,L::Ku);
    chem_state_replace(rng_, L::Ku, static_cast<size_t>(rhoP*nmols_), L::P);

    if (t0_active) {
      for (auto mi: *t0_active) {  if (chem_state_[mi]==L::Ku) {chem_state_[mi] = L::Kpp;}  }
    }
  }

protected:

  double match_events(mol_index_t mi1, mol_index_t mi2)
  {
    L const s1 = chem_state_[mi1],  s2 = chem_state_[mi2];

    if (     s1==L::Ku and s2==L::Ku) {
      return add_event(Ku_acts_on_Ku, mi1,mi2, 2*rate_c1);  // !Ku  + @Ku  --2*c1-->  Ku  + Kp
    }
    else if (s1==L::Ku and s2==L::Kp) {
      return add_event(Ku_acts_on_Kp, mi1,mi2, 1*rate_c1);  // !Ku  + @Kp  --1*c1-->  Ku  + Kpp
    }
    else if (s1==L::Kp and s2==L::Ku) {
      return add_event(Kp_acts_on_Ku, mi1,mi2, 2*rate_c2);  // !Kp  + @Ku  --2*c2-->  Kp  + Kp
    }
    else if (s1==L::Kp and s2==L::Kp) {
      return add_event(Kp_acts_on_Kp, mi1,mi2, 1*rate_c2);  // !Kp  + @Kp  --1*c2-->  Kp  + Kpp
    }
    else if (s1==L::Kpp and s2==L::Ku) {
      return add_event(Kpp_acts_on_Ku, mi1,mi2,2*rate_c3);  // !Kpp + @Ku  --2*c3-->  Kpp + Kp
    }
    else if (s1==L::Kpp and s2==L::Kp) {
      return add_event(Kpp_acts_on_Kp, mi1,mi2,1*rate_c3);  // !Kpp + @Kp  --1*c3-->  Kpp + Kpp
    }
    else if (s1==L::P and s2==L::Kpp) {
      return add_event(P_acts_on_Kpp, mi1,mi2, 2*rate_q0);  // !P   + @Kpp --2*c0-->  P   + Kp
    }
    else if (s1==L::P and s2==L::Kp) {
      return add_event(P_acts_on_Kp, mi1,mi2,  1*rate_q0);  // !P   + @Kp  --1*c0-->  P   + Ku
    }
    else { return 0.; }
  }

  static void Ku_acts_on_Ku (reagent_t& /*eKu*/, reagent_t& sKu) { *sKu.state()=L::Kp; }
  static void Ku_acts_on_Kp (reagent_t& /*eKu*/, reagent_t& sKp) { *sKp.state()=L::Kpp;}
  static void Kp_acts_on_Ku (reagent_t& /*eKp*/, reagent_t& sKu) { *sKu.state()=L::Kp; }
  static void Kp_acts_on_Kp (reagent_t& /*eKp*/, reagent_t& sKp) { *sKp.state()=L::Kpp;}
  static void Kpp_acts_on_Ku(reagent_t& /*eKpp*/,reagent_t& sKu) { *sKu.state()=L::Kp; }
  static void Kpp_acts_on_Kp(reagent_t& /*eKpp*/,reagent_t& sKp) { *sKp.state()=L::Kpp;}
  static void P_acts_on_Kpp (reagent_t& /*eP*/,  reagent_t& sKpp){*sKpp.state()=L::Kp; }
  static void P_acts_on_Kp  (reagent_t& /*eP*/,  reagent_t& sKp) { *sKp.state()=L::Ku; }

}; // ChemistrySystem2

#endif
