                                 __   __   __   ___  __
     _______________\\\    |    /  \ /  \ |__) |__  |__)    \\\_______________
    /               ///    |___ \__/ \__/ |    |___ |  \    ///               \
    |                                                                         |
    |  This software features the article "Sampling rare events in stochastic |
    |  reaction-diffusion systems within trajectory looping" by Zuk et al.,   |
    |  "Physical Review E" 98, 022401 (2018) DOI: 10.1103/PhysRevE.98.022401  |
    |                                                                         |
    |  Source code is distributed under the GNU GPL3 license.                 |
    |                                                                         |
    |  Homepage:  http://pmbm.ippt.pan.pl/software/looping  (permalink)       |
    \_________________________________________________________________________/


On the implementation
=====================

Trajectory looping can be implemented in two different but equivalent
ways.  One may first perform multiple joins of the base contacts
trajectory and in this manner obtain a looped contacts trajectory
of a desired length without considering chemical events.  The looped
contacts trajectory, consisting of the base contacts trajectory and
index maps resulting from optimal assignments, provides a complete
information necessary for a subsequent simulation of chemical
trajectory (this way has been presented in Fig. 1 in the article).

The second way consists in interweaving the simulation of chemical
events with joins of the base contacts trajectory.  This way can be
more convenient for some applications because what is expanding in
computer memory is the ultimate looped chemical trajectory, which is
of direct interest, and not the intermediary looped contacts trajectory,
whose desired length may be hard to guess.  When using the second way,
a stop condition can be defined depending on the chemical state of the
system or based on the already collected chemical event statistics.

This implementation of trajectory looping has been written according to
the second way.



Format of input trajectory
==========================

Looper reads in binary trajectories.  They are expected to have a header
consisting of:

* the number of molecules (4-byte signed integer),
* time step (double-precision real number),
* length of the edge of the cubical simulation box (positive value) or
  three (negative) values for lengths of non-cubical but cuboidal box
  (double-precision real number(s)).

The body of the trajectory file contains multiple frames. A single frame
consists of:

* current step index (4-byte signed integer),
* current time (double-precision real number),
* an ordered series of molecule coordinates.

Molecule coordinates comprise six values (six double-precision real numbers),
out of which three first describe absolute molecule position (_x_, _y_, _z_) and
the remaining three are only placeholders intended to be used in future to
describe molecule orientation. Currently, they are ignored (and thus can be
simply, e.g., three zeros).



Defining a new chemistry
========================

In order to change simulation parameters or define a new system of chemical
reactions, one has to modify and recompile the source code of Looper.

In file Chemistry.hpp there are two systems of reactions: `ChemistrySystem1`
and `ChemistrySystem2`, which correspond to the monostable and bistable
reaction systems, respectively, analyzed in the article.

To define and use a new system of reactions, one has to create a class that
inherits from the StochasticChemistry class.  The subclass should set initial
conditions in its constructor and should provide a `match_events(...)` method.
To provide a new system of reactions, please follow the conventions used to
implement `ChemistrySystem1` and `ChemistrySystem2` and associated enums for
molecule species.  The chemistry to be used in the simulator is pointed
in Settings.hpp with a type alias `chemistry_t`.



Various settings
================

Simulation duration
-------------------
Variable `kTimeEnd` in Settings.hpp.


Chemical rates and time scaling
-------------------------------
Variable `kLambda` in Settings.hpp.


Contact diameter
----------------
Variable `kSearchRadius` in Settings.hpp.


Skipping frames when writing output
-----------------------------------
Variable `kWriteToChemStateFileEverySteps` in Settings.hpp



Building
========

Compiler
--------
A decent and recent C++ compiler that supports C++14 standard is required.
The code has been confirmed to compile successfully under Linux with GNU
g++ 8.0 and Clang++ 6.0.


Compilation is managed by CMake.



Usage
=====

Invocation
----------
Invoke simply as:

    ./Looper trajectory chemprefix


Multiple input trajectories
---------------------------
More than one trajectory file can be read in and used.  Multiple files should
be listed consecutively before the last invocation argument, e.g., as:

    ./Looper traj1.bin traj2.bin traj3.bin chemprefix

The order of the use of multiple trajectories is random.  The trajectories
should have the same number of molecules and the same time step.


Multiple simulations: seed
--------------------------
If more than one simulation is to be performed (by setting appropriately
variable `kNumOfReplicas` in file Settings.hpp), each simulation is seeded
sequentially starting from the value assigned to `kSimulationSeedsBase` (file
Settings.hpp).  This value can be overridden by setting the environmental
variable `LOOPER_SEED_BASE` prior to Looper invocation.


Multiple simulations: multi-threading
-------------------------------------
Multiple simulations of chemical kinetics that use the same base chemical
trajectory can run in parallel.  During run-time, the software chooses the
number of threads based on both the `kNumOfReplicas` (described above) and
hardware capabilities.  If necessary, this behavior can be tweaked manually
using variable `kNumOfThreads` (in file Settings.hpp).



-- Last modified:  Wed Mar 28 23:44:17 CEST 2018

