/*
  Copyright (c) 2013 Uni Research AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "SimulatorOutput.hpp"

// we need complete definitions for these types
#include <opm/parser/eclipse/EclipseState/Schedule/TimeMap.hpp>
#include <opm/output/OutputWriter.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>

#include <numeric> // partial_sum

using namespace Opm;

SimulatorOutputBase::SimulatorOutputBase (
        const parameter::ParameterGroup& params,
        std::shared_ptr <const EclipseState> eclipseState,
        const Opm::PhaseUsage &phaseUsage,
        std::shared_ptr <const UnstructuredGrid> grid,
        std::shared_ptr <const SimulatorTimer> timer,
        std::shared_ptr <const SimulationDataContainer> state,
        std::shared_ptr <const WellState> wellState)

    // store all parameters passed into the object, making them curried
    // parameters to the writeOutput function.
    : timer_          (timer    )
    , reservoirState_ (state    )
    , wellState_      (wellState)

    // process parameters into a writer. we don't setup a new chain in
    // every timestep!
    , writer_ (std::move (OutputWriter::create (params, eclipseState, phaseUsage, grid)))
    // always start from the first timestep
    , next_ (0) {

    // write the static initialization files, even before simulation starts
    writer_->writeInit (*timer);
}

// default destructor is OK, just need to be defined
SimulatorOutputBase::~SimulatorOutputBase() { }

SimulatorOutputBase::operator std::function <void ()> () {
    // return (a pointer to) the writeOutput() function as an object
    // which can be passed to the event available from the simulator
    return std::bind (&SimulatorOutputBase::writeOutput, std::ref (*this));
}

void
SimulatorOutputBase::writeOutput () {
    const int this_time = timer_->simulationTimeElapsed ();

    // if the simulator signals for timesteps that aren't reporting
    // times, then ignore them
    if (next_ < timeMap_->size ()
        && timeMap_->getTimePassedUntil (next_) <= this_time) {
        // uh-oh, the simulator has skipped reporting timesteps that
        // occurred before this timestep (it doesn't honor the TSTEP setting)
        while (next_ < timeMap_->size ()
               && timeMap_->getTimePassedUntil (next_) < this_time) {
            ++next_;
        }

        // report this timestep if it matches
        if (next_ < timeMap_->size ()
            && timeMap_->getTimePassedUntil (next_) == this_time) {
            // make sure the simulator has spilled all necessary internal
            // state. notice that this calls *our* sync, which is overridden
            // in the template companion to call the simulator
            sync ();

            // relay the request to the handlers (setup in the constructor
            // from parameters)
            writer_->writeTimeStep (*timer_, *reservoirState_, *wellState_ , false);

            // advance to the next reporting time
            ++next_;
        }
    }
}

void
SimulatorOutputBase::sync () {
    // no-op in base class (overridden by simulator-specific template)
}
