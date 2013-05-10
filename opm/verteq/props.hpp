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
	 * @return Fluid object for the corresponding coarse grid. The caller
	 * has the responsibility to dispose off the object returned from here.
	 */
	static VertEqProps* create (const IncompPropertiesInterface& fineProps);
};

} // namespace Opm

#endif // OPM_VERTEQ_PROPS_HPP_INCLUDED
