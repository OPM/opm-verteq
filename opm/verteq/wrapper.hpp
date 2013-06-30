#ifndef OPM_VERTEQ_WRAPPER_HPP_INCLUDED
#define OPM_VERTEQ_WRAPPER_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#include <vector>

#ifndef OPM_VERTEQ_VISIBILITY_HPP_INCLUDED
#include <opm/verteq/visibility.hpp>
#endif /* OPM_VERTEQ_VISIBILITY_HPP_INCLUDED */

// forward declarations
struct UnstructuredGrid;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm {

namespace parameter { class ParameterGroup; }
class IncompPropertiesInterface;
class RockCompressibility;
class WellsManager;
class LinearSolverInterface;
class SimulatorTimer;
class TwophaseState;
class VertEq;
class WellsManager;
class WellState;
struct SimulatorReport;

/**
 * Wrapper that takes fine-scale 3D input, but uses an underlaying
 * simulator to do coarse-scale 2D simulation under the vertical
 * equilibrium assumption.
 *
 * The wrapper will setup the necessary translation layer for the
 * input, and provide methods to query a 3D translation of the output
 * (which may be expensive and is not done for every timestep).
 */
template <typename Simulator>
struct OPM_VERTEQ_PUBLIC VertEqWrapper {
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
	VertEqWrapper (
		const parameter::ParameterGroup& param,
		const UnstructuredGrid& grid,
		const IncompPropertiesInterface& props,
		const RockCompressibility* rock_comp_props,
		WellsManager& wells_manager,
		const std::vector<double>& src,
		const FlowBoundaryConditions* bcs,
		LinearSolverInterface& linsolver,
		const double* gravity);

	/**
	 * Clean up resources added in the wrapper.
	 */
	virtual ~VertEqWrapper ();

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

private:
	// vertical equilibrium model
	VertEq* ve;

	// list of translated wells
	WellsManager* wells_mgr;

	// underlaying simulator to use for 2D
	Simulator* sim;
};

} /* namespace Opm */

#endif /* OPM_VERTEQ_WRAPPER_HPP_INCLUDED */
