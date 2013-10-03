#ifndef OPM_VERTEQ_WRAPPER_HPP_INCLUDED
#define OPM_VERTEQ_WRAPPER_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#include <vector>
#include <memory> // unique_ptr

#ifndef OPM_VERTEQ_VISIBILITY_HPP_INCLUDED
#include <opm/verteq/visibility.hpp>
#endif /* OPM_VERTEQ_VISIBILITY_HPP_INCLUDED */

#ifndef OPM_VERTEQ_OPMFWD_HPP_INCLUDED
#include <opm/verteq/opmfwd.hpp>
#endif /* OPM_VERTEQ_OPMFWD_HPP_INCLUDED */

#ifndef OPM_VERTEQ_SIMULATOR_HPP_INCLUDED
#include <opm/verteq/simulator.hpp>
#endif /* OPM_VERTEQ_SIMULATOR_HPP_INCLUDED */

namespace Opm {

// forward declaration
class VertEq;

/**
 * Wrapper that takes fine-scale 3D input, but uses an underlaying
 * simulator to do coarse-scale 2D simulation under the vertical
 * equilibrium assumption.
 *
 * The wrapper will setup the necessary translation layer for the
 * input, and provide methods to query a 3D translation of the output
 * (which may be expensive and is not done for every timestep).
 *
 * You are not supposed to use the base class directly; rather use
 * VertEqWrapper<T>.
 */
struct OPM_VERTEQ_PUBLIC VertEqWrapperBase {
protected:
	/**
	 * Initialise from parameters and objects to observe.
	 *
	 * @param underlaying     Type to use for 2D simulations. This pointer
	 *                        is adopted, i.e. the wrapper takes ownership
	 *                        of it. Use with newly created SimulatorInstance.
	 * @param param           Parameters for the underlaying simulator class
	 * @param grid            Fine-scale grid data structure
	 * @param props           Fluid and rock properties
	 * @param rock_comp_props If non-null, rock compressibility properties
	 * @param well_manager    Well manager, may manage no (null) wells
	 * @param src             Source terms
	 * @param bcs             Boundary conditions, treat as all noflow if null
	 * @param linsolver       Linear solver
	 * @param gravity         If non-null, gravity vector
	*/
	VertEqWrapperBase (
		std::unique_ptr <Simulator> underlaying,
		const parameter::ParameterGroup& param,
		const UnstructuredGrid& grid,
		const IncompPropertiesInterface& props,
		const RockCompressibility* rock_comp_props,
		WellsManager& wells_manager,
		const std::vector<double>& src,
		const FlowBoundaryConditions* bcs,
		LinearSolverInterface& linsolver,
		const double* gravity);
public:
	/**
	 * Clean up resources added in the wrapper.
	 */
	virtual ~VertEqWrapperBase ();

	/**
	 * Run the simulation.
	 *
	 * This will run succesive timesteps until timer.done() is true. It will
	 * modify the reservoir and well states.
	 *
	 * @param[in,out] timer       Governs the requested reporting timesteps
	 * @param[in,out] state       State of reservoir: pressure, fluxes
	 * @param[in,out] well_state  State of wells: bhp, perforation rates
	 * @return                    Simulation report, with timing data
	 */
  virtual SimulatorReport run(
		SimulatorTimer& timer,
		TwophaseState& state,
		WellState& well_state);

	/**
	 * Event that is signaled every time the simulator has completed a
	 * timestep.
	 *
	 * @return Reference to the event object where callbacks can be registered.
	 *
	 * @see Opm::SimulatorIncompTwophase::timestep_completed
	 */
	Event& timestep_completed ();

	/**
	 * Notify the simulator that a callback has an interest in reading
	 * for reporting purposes the contents of the state argument that
	 * was passed to the run() method. The simulator will then flush
	 * any internal state which is currently not reflected in it.
	 *
	 * @see Opm::SimulatorIncompTwophase::sync
	 */
	void sync ();

private:
	// underlaying simulator to use for 2D
	std::unique_ptr <Simulator> sim;

	// vertical equilibrium model
	VertEq* ve;

	// list of translated wells
	WellsManager* wells_mgr;

	// list of handlers for timestep notifications
	std::unique_ptr <EventSource> timestep_callbacks;

	// current state of any simulation we are doing
	TwophaseState* fineState;
	TwophaseState* coarseState;

	// flag that determines whether we have synced or not
	bool syncDone;
	void resetSyncFlag ();
};

/**
 * Wrap a certain type of simulator to use for in an vertically upscaled
 * model.
 *
 * This is its own class just because you cannot in C++ call a templated
 * constructor without the template as any of the parameters. (The
 * template constructor sets up the adapter class). By using an adapter,
 * we don't have to instantiate the wrapper for every simulator class.
 *
 * @tparam SimulatorType Type of the underlaying simulator to use for
 *                       2D simulations, e.g. SimulatorIncompTwophase.
 */
template <typename SimulatorType>
struct OPM_VERTEQ_PUBLIC VertEqWrapper : public VertEqWrapperBase {
	/**
	 * Initialise from parameters and objects to observe.
	 *
	 * @param param           Parameters for the underlaying simulator class
	 * @param grid            Fine-scale grid data structure
	 * @param props           Fluid and rock properties
	 * @param rock_comp_props If non-null, rock compressibility properties
	 * @param well_manager    Well manager, may manage no (null) wells
	 * @param src             Source terms
	 * @param bcs             Boundary conditions, treat as all noflow if null
	 * @param linsolver       Linear solver
	 * @param gravity         If non-null, gravity vector
	 */
	VertEqWrapper (const parameter::ParameterGroup& param,
	               const UnstructuredGrid& grid,
	               const IncompPropertiesInterface& props,
	               const RockCompressibility* rock_comp_props,
	               WellsManager& wells_manager,
	               const std::vector<double>& src,
	               const FlowBoundaryConditions* bcs,
	               LinearSolverInterface& linsolver,
	               const double* gravity)
		: VertEqWrapperBase (std::unique_ptr <Simulator> (
		                       new SimulatorAdapter <SimulatorType> ()),
		                     param, grid, props, rock_comp_props, wells_manager,
		                     src, bcs, linsolver, gravity) {
	}
};

} /* namespace Opm */

#endif /* OPM_VERTEQ_WRAPPER_HPP_INCLUDED */
