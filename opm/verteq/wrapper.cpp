#include <opm/verteq/wrapper.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
using namespace Opm;

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
	: sim (0) {

	// pass arguments to the underlaying simulator
	sim = new Simulator (param,
	                     grid,
	                     props,
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
}

template <typename Simulator>
SimulatorReport
VertEqWrapper <Simulator>::run(
		SimulatorTimer& timer,
		TwophaseState& state,
    WellState& well_state) {

	// forward the call to the underlaying simulator
	return sim->run(timer, state, well_state);
}

// supported wrappings; we need to list here all the possible
// instantiations that should be done and put in the library
// (and no other will be possible because they are missing the
// code which is provided in this compilation unit).
#include <opm/core/simulator/SimulatorIncompTwophase.hpp>
template class VertEqWrapper <SimulatorIncompTwophase>;
