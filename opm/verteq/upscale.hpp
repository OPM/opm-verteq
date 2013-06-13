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

protected:
	const TopSurf& ts;
};

} // namespace Opm

#endif // OPM_VERTEQ_UPSCALE_HPP_INCLUDED
