#ifndef OPM_VERTEQ_TOPSURF_HPP_INCLUDED
#define OPM_VERTEQ_TOPSURF_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef BOOST_RANGE_ITERATOR_RANGE_INCLUDED
#include <boost/range/iterator_range.hpp>
#endif

#ifndef OPM_GRID_HEADER_INCLUDED
#include <opm/core/grid.h>
#endif

namespace Opm {

/**
 * Two-dimensional top surface of a full, three-dimensional grid.
 *
 * This grid is set up such that each cell is an upscaling of all the cells
 * in each column in the full grid. It also contains a mean to map results
 * from this grid back to the full grid.
 *
 * The full grid is also referred to as the fine grid, and this grid as the
 * coarse grid, or upscaled, grid.
 *
 * Note: Do NOT call destroy_grid () when done with this structure; it will
 * only clean up half of it. Wrap it is a smart pointer that calls the
 * destructor.
 */
struct TopSurf : public UnstructuredGrid {
	virtual ~TopSurf ();

	/**
	 * Indices of the columns' underlaying cells in the full grid.
	 *
	 * Consecutive indices from the _fine_ grid, not this one, for each of the
	 * columns, i.e. cells in _this_ grid, as one flat list.
	 *
	 * The values are sorted in z-order, starting from the top and moving
	 * downwards to the bottom. This is useful because you can keep a running
	 * counter for the depth, filling items as you go. (For this to be really
	 * useful, the original grid should be reordered so that cells in the
	 * z-direction are closer).
	 *
	 * Use this field together with the col_cellpos to iterate through a column
	 * in the fine grid.
	 *
	 * @see TopSurf::column, TopSurf::col_cellpos
	 */
	int* col_cells;

	/**
	 * Number of cells in the columns preceeding each one.
	 *
	 * For each column c, the number col_cellpos[c] is the number of cells in
	 * the _full_ grid that belongs to the columns 0..(c-1).
	 *
	 * This arrangement means that col_cellpos[c] is the index into col_cells
	 * of the first fine cell, whereas col_cellpos[c+1] is the index into
	 * col_cells of the last fine cell.
	 *
	 * @see TopSurf::column, TopSurf::col_cellpos
	 */
	int* col_cellpos;

	/**
	 * Create an upscaled grid based on a full, three-dimensional grid.
	 *
	 * @param fine Grid that should be upscaled.
	 *
	 * This must be a three-dimensional, Cartesian grid which has an active
	 * cluster which is without holes and which is convex. The UnstructuredGrid
	 * structure is used because it is the lingua franca of grids in the
	 * simulator, not because this method will handle every possible grid.
	 *
	 * This pointer is NOT adopted. The caller must still dispose of the grid.
	 *
	 * @return Upscaled, fine grid.
	 *
	 * The caller have the responsibility of disposing this grid; no other
	 * references will initially exist.
	 */
	static TopSurf* create (const UnstructuredGrid& fine);

	/**
	 * Iterator for column members
	 *
	 * Defines a input iterator that let one walk through all the cells in
	 * the fine grid, which belongs to a particular column, i.e. cell in the
	 * coarse grid.
	 *
	 * Each item pointed to by this iterator is an index into the full grid.
	 * The iterator is only valid for one column.
	 *
	 * @param ndx2d Index into the two-dimensional grid.
	 *
	 * Each of the elements in the two-dimensional grid identifies a column
	 * in the original, full grid.
	 *
	 * @return Iterator through three-dimensional elements in the column.
	 *
	 * The iterator supports the idiomatic methods begin() and end() which
	 * can be used to get to the starting and ending point of the iteration,
	 * respectively.
	 *
	 * @example
	 * @code{.cpp}
	 * BOOST_FOREACH (int ndx_3d, topSurf->column (ndx_2d)) {
	 *   sum += poro [ndx_3d];
	 * }
	 * @endcode
	 */
	boost::iterator_range <const int*> column (int ndx_2d);

private:
	/**
	 * @brief You are not meant to construct these yourself; use create ().
	 */
	TopSurf ();
};

} // namespace Opm

#endif // OPM_VERTEQ_TOPSURF_HPP_INCLUDED
