#ifndef OPM_VERTEQ_VERTEQ_HPP_INCLUDED
#define OPM_VERTEQ_VERTEQ_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_VERTEQ_VISIBILITY_HPP_INCLUDED
#include <opm/verteq/visibility.hpp>
#endif /* OPM_VERTEQ_VISIBILITY_HPP_INCLUDED */

// forward declaration
struct UnstructuredGrid;

namespace Opm {

namespace parameter {
struct ParameterGroup;
} // namespace parameter

/**
 * @brief Vertical Equilibrium Upscaling
 *
 * This class acts like a wrapper around a three-dimensional setup of
 * grid, fluid and rock properties and provide a view of a two-
 * dimensional upscaling based on the vertical equilibrium model.
 *
 * The objects provided by the upscaling may be passed into the
 * simulator as if read from a deck, and the results translated back
 * into a full model afterwards.
 *
 * Notice that the following constraints apply on the grid:
 *
 * * It must be structured, i.e. every element is (i,j,k)-addressable
 * * No degenerate faces (on the top), i.e. no crossing grid axes
 * * It must have no horizontal faults
 */
class OPM_VERTEQ_PUBLIC VertEq
{
protected:
	// you are not supposed to call the constructor yourself; instead use
	// the static create function provided
	VertEq () {}
public:
	/**
	 * @brief Pseudo-constructor of VertEqUpscaling objects
	 *
	 * You own the object that is returned, and is responsible for
	 * calling delete on the pointer.
	 *
	 * @param title Name of the case, gotten from getTITLE().name(); this
	 *              may be used to set grid-specific properties.
	 * @param args Parameters
	 * @param fullGrid Grid obtained elsewhere. This object is not
	 *        adopted, but is assumed to be live over the lifetime
	 *        of the upscaling.
	 *
	 * @return A new instance implementing this interface.
	 */
	static VertEq* create (const std::string& title,
												 const Opm::parameter::ParameterGroup& args,
												 const UnstructuredGrid& fullGrid);

	// virtual destructor, actual functionality relayed to real impl.
	virtual ~VertEq () {}

	/**
	 * @brief Accessor method for the upscaled grid
	 *
	 * This method is inexpensive; the grid is not constructed upon
	 * every invocation.
	 *
	 * @return Grid that may be passed to other components in the
	 *         simulator. You do NOT own this object!
	 *
	 */
	virtual const UnstructuredGrid& grid () = 0;
};

} // namespace Opm

#endif /* OPM_VERTEQ_VERTEQ_HPP_INCLUDED */
