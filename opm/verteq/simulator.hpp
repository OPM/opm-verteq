#ifndef OPM_VERTEQ_SIMULATOR_HPP_INCLUDED
#define OPM_VERTEQ_SIMULATOR_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#include <memory> // unique_ptr
#include <vector>

#ifndef OPM_VERTEQ_VISIBILITY_HPP_INCLUDED
#include <opm/verteq/visibility.hpp>
#endif /* OPM_VERTEQ_VISIBILITY_HPP_INCLUDED */

#ifndef OPM_SIMULATORREPORT_HEADER_INCLUDED
#include <opm/core/simulator/SimulatorReport.hpp>
#endif /* OPM_SIMULATORREPORT_HEADER_INCLUDED */

#ifndef OPM_VERTEQ_OPMFWD_HPP_INCLUDED
#include <opm/verteq/opmfwd.hpp>
#endif /* OPM_VERTEQ_OPMFWD_HPP_INCLUDED */

namespace Opm {

/**
 * Interface to which the underlaying two-phase simulator must adhere.
 *
 * Notice that it must support delayed construction through the init
 * method which implicates that implementations of this will probably
 * be some kind of smart pointer wrapper.
 */
struct OPM_VERTEQ_PUBLIC Simulator {
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
	virtual void init (
		const parameter::ParameterGroup& param,
		const UnstructuredGrid& grid,
		const IncompPropertiesInterface& props,
		const RockCompressibility* rock_comp_props,
		WellsManager& wells_manager,
		const std::vector<double>& src,
		const FlowBoundaryConditions* bcs,
		LinearSolverInterface& linsolver,
		const double* gravity) = 0;

	/**
	 * Clean up resources added in the wrapper.
	 */
    virtual ~Simulator () {}

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
		WellState& well_state) = 0;

	/**
	 * Event that is signaled every time the simulator has completed a
	 * timestep.
	 *
	 * @return Reference to the event object where callbacks can be registered.
	 *
	 * @see Opm::SimulatorIncompTwophase::timestep_completed
	 */
	virtual Event& timestep_completed () = 0;

	/**
	 * Notify the simulator that a callback has an interest in reading
	 * for reporting purposes the contents of the state argument that
	 * was passed to the run() method. The simulator will then flush
	 * any internal state which is currently not reflected in it.
	 *
	 * @see Opm::SimulatorIncompTwophase::sync
	 */
	virtual void sync () = 0;
};

/**
 * Instance of a concrete simulator that adhere to the interface. This
 * adapter class allows us to bring into the simulator type class any
 * class that has the necessary methods but not necessarily shares a
 * common base class.
 *
 * @tparam T Type of the actual simulator, e.g. SimulatorIncompTwophase
 */
template <typename T>
struct OPM_VERTEQ_PUBLIC SimulatorAdapter : public Simulator {
    /**
     * Activate the reference by creating a new instance of the
     * underlaying simulator with the given parameters and then
     * point to it.
     */
	virtual void init (
		const parameter::ParameterGroup& param,
		const UnstructuredGrid& grid,
		const IncompPropertiesInterface& props,
		const RockCompressibility* rock_comp_props,
		WellsManager& wells_manager,
		const std::vector<double>& src,
		const FlowBoundaryConditions* bcs,
		LinearSolverInterface& linsolver,
		const double* gravity) {
        t_ = std::unique_ptr <T> (new T (param, grid, props, rock_comp_props,
                    wells_manager, src, bcs, linsolver,
                    gravity));
    }

    // forward all method to the underlaying instance
    virtual SimulatorReport run(
		SimulatorTimer& timer,
		TwophaseState& state,
		WellState& well_state) {
        return t_->run (timer, state, well_state);
    }
    virtual Event& timestep_completed () {
        return t_->timestep_completed ();
    }
    virtual void sync () {
        t_->sync ();
    }

private:
    std::unique_ptr <T> t_;
};

} /* namespace Opm */

#endif /* OPM_VERTEQ_SIMULATOR_HPP_INCLUDED */
