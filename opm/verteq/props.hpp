#ifndef OPM_VERTEQ_PROPS_HPP_INCLUDED
#define OPM_VERTEQ_PROPS_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_VERTEQ_VISIBILITY_HPP_INCLUDED
#include <opm/verteq/visibility.hpp>
#endif /* OPM_VERTEQ_VISIBILITY_HPP_INCLUDED */

#ifndef OPM_INCOMPPROPERTIESINERFACE_HEADER_INCLUDED
#include <opm/core/props/IncompPropertiesInterface.hpp>
#endif /* OPM_INCOMPPROPERTIESINERFACE_HEADER_INCLUDED */

namespace Opm {

// forward declarations
struct TopSurf;

struct VertEqProps : public IncompPropertiesInterface {
	/**
	 * Create an upscaled version of fluid and rock properties.
	 *
	 * @param fineProps Fluid and rock properties for the fine grid.
	 * @param topSurf Grid for which the properties should be upscaled.
	 * @param gravity Gravity vector (three-dimensional); must contain
	 *                three elements, whereas the last is for depth.
	 *                Usually this is {0., 0., Opm::unit::gravity}.
	 * @return Fluid object for the corresponding coarse grid. The caller
	 * has the responsibility to dispose off the object returned from here.
	 */
	static VertEqProps* create (const IncompPropertiesInterface& fineProps,
	                            const TopSurf& topSurf,
	                            const double* gravity);

	/**
	 * Update residual saturation of CO2 through-out the domain.
	 *
	 * When the plume move forward, it leaves behind some residual
	 * CO2 which wasn't there before. This must be accounted for when
	 * calculating the position of the interface from the upscaled
	 * saturation.
	 *
	 * @param sat Saturation for the phases in the coarse domain. It
	 *            consists of one record for each column in the top
	 *            surface grid, where each record has one saturation
	 *            for each phase. The ordering of the phases is the
	 *            same as in the properties interface.
	 */
	virtual void upd_res_sat (const double* sat) = 0;

	/**
	 * Upscale pressure from fine-scale to coarse-scale.
	 *
	 * @param[in] coarseSaturation
	 *	Saturation for each phase, in each column. This is used to initialize
	 *	the brine-co2 phase contact properly.
	 *
	 * @param[in] finePressure
	 *	Pressure for each block in the fine grid. There is only one value
	 *	for pressure, the pressure for hydrostatic equilibrium in each block.
	 *
	 * @param[out] coarsePressure
	 *	Pressure for each column in the top surface grid. There is only one
	 *	value for pressure, weighted from each of the blocks in the column.
	 *
	 *	The space for the data must have been allocated by the caller.
	 */
	virtual void upscale_pressure (const double* coarseSaturation,
	                               const double* finePressure,
	                               double* coarsePressure) = 0;

	/**
	 * Upscale saturation from fine-scale to coarse-scale.
	 *
	 * @param[in] fineSaturation
	 *	Saturation for each phase, and for each block in the fine grid. The
	 *	data for each block is kept together, i.e. the phase index varies
	 *	most quickly.
	 *
	 * @param[out] coarseSaturation
	 *	Saturation for each phase, and for each column in the coarse grid. The
	 *	data for each column is kept together, i.e. the phase index varies
	 *	most quickly.
	 *
	 *	The space for the data must have been allocated by the caller.
	 */
	virtual void upscale_saturation (const double* fineSaturation,
	                                 double* coarseSaturation) = 0;

	/**
	 * Downscale to corresponding 3D fine-scale saturations from 2D
	 * coarse-scale saturations.
	 *
	 * @param[in] coarseSaturation
	 *	Saturation for each phase, and for each column in the coarse grid. The
	 *	data for each column is kept together, i.e. the phase index varies
	 *	most quickly.
	 *
	 * @param[out] fineSaturation
	 *	Saturation for each phase, and for each block in the fine grid. The
	 *	data for each block is kept together, i.e. the phase index varies
	 *	most quickly.
	 *
	 * @note The space for the data must have been allocated by the caller.
	 */
	virtual void downscale_saturation (const double* coarseSaturation,
	                                   double* fineSaturation) = 0;

	/**
	 * Downscale to corresponding 3D fine-scale pressure from 2D coarse-scale
	 * pressure.
	 *
	 * @param[in] coarseSaturation
	 *	Saturation for each phase, in each column. This is used to determine
	 *	the brine-co2 phase contact properly.
	 *
	 * @param[in] coarsePressure
	 *	Pressure of the CO2 phase at the top of each column.
	 *
	 * @param[out] finePressure
	 *	Pressure in each block in the fine-scale grid. The order of the cells
	 *	is given by the fine_col/col_cells members of the TopSurf object which
	 *	was passed the constructor.
	 *
	 * @note The space for the data must have been allocated by the caller.
	 */
	virtual void downscale_pressure (const double* coarseSaturation,
	                                 const double* coarsePressure,
	                                 double* finePressure) = 0;
};

} // namespace Opm

#endif // OPM_VERTEQ_PROPS_HPP_INCLUDED
