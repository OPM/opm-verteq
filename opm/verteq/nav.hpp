#ifndef OPM_VERTEQ_NAV_HPP_INCLUDED
#define OPM_VERTEQ_NAV_HPP_INCLUDED

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

#endif // OPM_VERTEQ_NAV_HPP_INCLUDED
