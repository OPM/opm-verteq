#ifndef OPM_VERTEQ_UPSCALE_HPP_INCLUDED
#define OPM_VERTEQ_UPSCALE_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_VERTEQ_VISIBILITY_HPP_INCLUDED
#include <opm/verteq/visibility.hpp>
#endif /* OPM_VERTEQ_VISIBILITY_HPP_INCLUDED */

#ifndef OPM_VERTEQ_TOPSURF_HPP_INCLUDED
#include <opm/verteq/topsurf.hpp>
#endif /* OPM_VERTEQ_TOPSURF_HPP_INCLUDED */

namespace Opm {

/**
 * Extension of the top surface that does integration across columns.
 * The extension is done by aggregation since these methods really are
 * orthogonal to how the grid is created.
 */
struct VertEqUpscaler {
	/**
	 * Create a layer on top of a top surface that can upscale properties
	 * in columns.
	 *
	 * @param topSurf Grid that provides the weights in the integrations.
	 */
	VertEqUpscaler (const TopSurf& topSurf)
	  : ts (topSurf) {
	}

	/**
	 * Retrieve a property from the fine grid for a certain column of blocks.
	 *
	 * The property is stored in records of 'stride' length, 'offset' elements
	 * from the start. (This is used to index into a permeability tensor).
	 *
	 * @param col Index of the column to retrieve. This must be a valid cell
	 *            number in the top surface.
	 *
	 * @param buf Array that will receive the data. It is assumed that this
	 *            is preallocated, and is large enough to hold the data for
	 *            this column.
	 *
	 * @param data Array that contains the property for the *entire* fine grid.
	 *             This is typically something that was retrieved from a fluid
	 *             object associated with the fine grid.
	 *
	 * @param stride Number of values (not bytes!) between each entry for the
	 *               property. Use one for simple arrays.
	 *
	 * @param offset Number of values (not bytes!) before the first entry for
	 *               this property. Use zero for simple arrays.
	 */
	void gather (int col, double* buf, const double* data, int stride, int offset) const;

	/**
	 * Depth fraction weighted by an expression specified as piecewise values
	 * downwards a column from the top.
	 *
	 * Returned for each block is the integral from the top down to and
	 * including each block, divided by the total height of the column. The
	 * last value will thus contain the average for the entire column.
	 *
	 * @param col Index of the column for which the expression should be
	 *            integrated.
	 *
	 * @param val Integrand values for each block in the column.
	 *
	 * @param res Array that will receive the result. The caller must
	 *            preallocate storage and pass it as input here. A good idea
	 *            is to use the view of a column from an rlw_double.
	 */
	void wgt_dpt (int col, const double* val, double* res) const;

protected:
	const TopSurf& ts;
};

} // namespace Opm

#endif // OPM_VERTEQ_UPSCALE_HPP_INCLUDED
