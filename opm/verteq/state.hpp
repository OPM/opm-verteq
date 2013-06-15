#ifndef OPM_VERTEQ_STATE_HPP_INCLUDED
#define OPM_VERTEQ_STATE_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_VERTEQ_VISIBILITY_HPP_INCLUDED
#include <opm/verteq/visibility.hpp>
#endif /* OPM_VERTEQ_VISIBILITY_HPP_INCLUDED */

#ifndef OPM_TWOPHASESTATE_HEADER_INCLUDED
#include <opm/core/simulator/TwophaseState.hpp>
#endif /* OPM_TWOPHASESTATE_HEADER_INCLUDED */

#ifndef OPM_VERTEQ_VERTEQ_HPP_INCLUDED
#include <opm/verteq/verteq.hpp>
#endif /* OPM_VERTEQ_VERTEQ_HPP_INCLUDED */

namespace Opm {

/**
 * State object that knows which model that created it and which
 * can communicate changes back to it.
 *
 * Since the state is not part of the construction of the simulator
 * object, it is not available when we need to create the VE model,
 * but must be supplied later.
 *
 * Use this object as a replacement for the fine-scale state.
 */
struct VertEqState : public TwophaseState {
	VertEq& ve;

	VertEqState (VertEq& model,
	             const TwophaseState& fineScale)
		// keep a pointer back to the model, so we know who to call
		: ve (model) {

		// let the VE model do the initialization, it has all the
		// necessary information to do the upscaling
		ve.upscale (fineScale, *this);
	}

	/**
	 * Push changes that are done to this object back into the VE
	 * model when it has been updated (a timestep has completed)
	 *
	 * (The reason this method isn't const is because it logically
	 * changes the underlaying bundled VE model)
	 *
	 * @see SimulatorIncompTwophase::connect_timestep
	 */
	void notify () {
		ve.notify (*this);
	}
};

} /* namespace Opm */

#endif /* OPM_VERTEQ_STATE_HPP_INCLUDED */
