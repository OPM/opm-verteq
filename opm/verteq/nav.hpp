#ifndef OPM_VERTEQ_NAV_HPP_INCLUDED
#define OPM_VERTEQ_NAV_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_GRID_HEADER_INCLUDED
#include <opm/core/grid.h>
#endif

/**
 * There are three types of indices used in this module:
 *
 * (1) Cartesian coordinate
 * (2) Cartesian index
 * (3) Global identity
 *
 * The Cartesian coordinate is an (i,j)-tuple into the cornerpoint grid
 * structure.
 *
 * The Cartesian index, is a flat integer which has is
 * determined solely by the structure of the grid, regardless of whether
 * there are any active elements. It can be calculated from the coordinates
 * if the extent of the grid is known. We use this to efficiently store
 * data for each possible position in the grid without using a predefined
 * size or dynamic structure for each column.
 *
 * The global identity is the index which is assigned to it in the grid
 * structure, after inactive elements are discarded. This is the 'identity'
 * of the cell in the rest of the simulator. Cells that aren't active are
 * not assigned an identity.
 *
 * The value types defined here provide a way to address location in
 * the grid in a type-safe manner to let the compiler help us keep track
 * of the real meaning of indices.
 *
 * The navigation classes provide a way to define an enumeration of the
 * grid without resorting to inline integer arithmetic inside the other
 * functions. An optimizing compiler should be able to generate equally
 * fast code as hand-coded index calculations, when using these classes.
 */

/**
 * Index tuple in two-dimensional cornerpoint grid
 *
 * This structure represents the carrier of Cartesian coordinates. They
 * should be thought of as an integral type, along the lines of complex
 * numbers.
 */
struct Coord2D {
	const int i;
	const int j;

	Coord2D (int i_, int j_) : i (i_), j (j_) {	}
};

/// Index tuple in three-dimensional cornerpoint grid.
struct Coord3D : public Coord2D {
	const int k;

	Coord3D (int i_, int j_, int k_)
		: Coord2D (i_, j_)
		, k (k_) {
	}
};

/// Type-safe enumeration of axis directions.
struct Dir {
	/// Towards the end of the axis with lesser numbers.
	static const Dir DEC; // = 0

	/// Towards the end of the axis with greater numbers.
	static const Dir INC; // = 1

	/// Number of possible directions
	static const int COUNT = 2;

	/// Integer representation suitable for indexing in array
	const int val;

	Dir (const Dir& rhs) : val (rhs.val) {}

protected:
	/// Private constructor to avoid initialization outside domain
	Dir (int i) : val (i) { }
};

/// Type-safe enumeration of axis dimensions
struct Dim2D {
	// two spatial directions
	static const Dim2D X; // = 0
	static const Dim2D Y; // = 1

	// number of dimensions
	static const int COUNT = 2;

	const int val;

	Dim2D (const Dim2D& rhs) : val (rhs.val) { }

protected:
	Dim2D (int i) : val (i) { }
};

/// Value type that addresses sides in a two-dimensional grid cell
struct Side2D {
	const Dim2D dim;
	const Dir dir;

	Side2D (Dim2D dim_, Dir dir_) : dim (dim_), dir (dir_) { }
	Side2D (const Side2D& rhs) : dim (rhs.dim), dir (rhs.dir) {}
};

/// Value type that addresses corners in a two-dimensional grid cell
struct Corn2D {
	const Dir i;
	const Dir j;
	Corn2D (Dir i_, Dir j_) : i (i_), j (j_) { }
	Corn2D (const Corn2D& rhs) : i (rhs.i), j (rhs.j) { }
};

/**
 * Navigate a Cartesian grid in a structured way so that clearly defined
 * mapping between the enumeration index and the coordinate.
 */
struct Cart2D {
	// number of cells in each direction
	const int ni;
	const int nj;

	// initialize from the size of the grid
	Cart2D (int ni_, int nj_) : ni (ni_), nj (nj_) { }

	// use these two value types for structured coordinates and
	// flattened indices, respectively
	typedef Coord2D coord_t;
	typedef int elem_t;

	/// Number of (possible) elements in the grid
	int num_elems () const {
		return ni * nj;
	}

	/// Cartesian (flattened) index for a coordinate
	elem_t cart_ndx (const coord_t& coord) const {
		return coord.j * ni + coord.i;
	}

	/// Cartesian coordinate for a (flattened) index
	coord_t coord (const elem_t& cart_ndx) const {
		const div_t strip = div (cart_ndx, ni);
		const int i = strip.rem;
		const int j = strip.quot;
		return coord_t (i, j);
	}
};

/**
 * Navigate a three-dimensional grid.
 *
 * In this module, we only need this to get the structured index of
 * each three-dimensional element, in order to know into which column
 * we should assign this element. However, we keep the design in the
 * same manner as the two-dimensional case.
 */
struct Cart3D {
	// number of cells in each direction
	const int ni;
	const int nj;
	const int nk;

	/// Initialize POD from an existing (3D) grid
	Cart3D (const UnstructuredGrid& g)
		: ni (g.cartdims [0])
		, nj (g.cartdims [1])
		, nk (g.cartdims [2]) { }

	/// Project grid into a surface
	Cart2D project () const {
		return Cart2D (ni, nj);
	}

	// use these two value types for structured coordinates and
	// flattened indices, respectively
	typedef Coord3D coord_t;
	typedef int elem_t;

	/// Deconstruct Cartesian index into coordinates
	coord_t coord (const elem_t& cart_ndx) const {
		// the i-index moves fastest, as this is Fortran-indexing
		const div_t strip = div (cart_ndx, ni);
		const int i = strip.rem;
		const div_t plane = div (strip.quot, nj);
		const int j = plane.rem;
		const int k = plane.quot;
		return coord_t (i, j, k);
	}
};

#endif // OPM_VERTEQ_NAV_HPP_INCLUDED
