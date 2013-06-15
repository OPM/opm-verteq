#include <opm/verteq/wrapper.hpp>
#include <string>
#include <opm/verteq/verteq.hpp>
#include <opm/verteq/state.hpp>
#include <opm/core/simulator/SimulatorIncompTwophase.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
using namespace std;

namespace Opm {

template <typename Simulator>
VertEqWrapper <Simulator>::VertEqWrapper (
		const parameter::ParameterGroup& param,
		const UnstructuredGrid& grid,
		const IncompPropertiesInterface& props,
		const RockCompressibility* rock_comp_props,
		WellsManager& wells_manager,
		const std::vector<double>& src,
		const FlowBoundaryConditions* bcs,
		LinearSolverInterface& linsolver,
		const double* gravity)

	// initialize pointer to null, so it is deleted properly
	// if the sim. ctor throws
	: ve (0)
	, sim (0) {

	// VE model that is injected in between the fine-scale
	// model that is sent to us, and the simulator
	ve = VertEq::create ("",
	                     param,
	                     grid,
	                     props);

	// pass arguments to the underlaying simulator
	sim = new Simulator (param,
	                     ve->grid (),
	                     ve->props (),
	                     rock_comp_props,
	                     wells_manager,
	                     src,
	                     bcs,
	                     linsolver,
	                     gravity);
}

template <typename Simulator>
VertEqWrapper <Simulator>::~VertEqWrapper () {
	delete sim;
	delete ve;
}

// need to template instantiate this method in order to use the
// template method within it (connect_timestep); this is not a
// disaster; different simulators will probably need different
// specializations
template <> SimulatorReport
VertEqWrapper <SimulatorIncompTwophase>::run(
		SimulatorTimer& timer,
		TwophaseState& state,
		WellState& well_state) {

	// the state that is passed is for fine-scale; we need to create
	// a coarse state that matches the grid, for the simulator
	VertEqState upscaled_state (*ve, state);

	// make the state "active", so that it push its changes to the
	// ve model whenever an update is completed and its state is stable
	sim->connect_timestep <VertEqState, &VertEqState::notify> (upscaled_state);

	// forward the call to the underlaying simulator
	return sim->run(timer, upscaled_state, well_state);
}

// supported wrappings; we need to list here all the possible
// instantiations that should be done and put in the library
// (and no other will be possible because they are missing the
// code which is provided in this compilation unit).
template class VertEqWrapper <SimulatorIncompTwophase>;

} /* namespace Opm */
