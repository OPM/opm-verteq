// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#include <opm/verteq/nav.hpp>
#include <opm/verteq/topsurf.hpp>
#include <opm/core/utility/exc.hpp>
#include <boost/range/iterator_range.hpp>
#include <algorithm> // min, max
#include <climits> // INT_MIN, INT_MAX
#include <cstdlib> // div
#include <map>
#include <memory> // auto_ptr
#include <utility> // pair

using namespace boost;
using namespace Opm;
using namespace std;

/**
 * Generic run-length iterator.
 *
 * @param ndx Index of the range to find, e.g. a column.
 *
 * @param pos Array containing accumulated counts for each range,
 * e.g. col_collpos.
 *
 * @param values Array containing the values for every range,
 * concatenated into one huge array, e.g. col_cells.
 *
 * @return Iterator that can be used to iterate through all values
 * only in the specified range, e.g. cells in a column.
 */
template <typename T> iterator_range <const T*>
run_len_iter (const int ndx, const int* const& pos, const T* const& values) {
	// skip all the values that belongs to ranges before us,
	// then we get to the start of our own range
	const T* begin_addr = &values [pos [ndx]];

	// stop when we arrive at the start of the next range
	const T* end_addr = &values [pos [ndx + 1]];

	// return an iterator over this
	return iterator_range <const T *> (begin_addr, end_addr);
}

iterator_range <const int*>
TopSurf::column (int ndx_2d) {
	return run_len_iter <const int> (ndx_2d, this->col_cellpos, this->col_cells);
}

/**
 * @brief Process to extract the top surface from a structured grid.
 *
 * This object encapsulates a procedure with variables shared amongst
 * several sub-procedures (like in Pascal). These objects are not
 * supposed to linger on afterwards.
 */
struct TopSurfBuilder {
	// source grid from which we get the input data
	const UnstructuredGrid& fine_grid;

	// target grid we are constructing
	TopSurf& ts;

	// number of grid dimensions in each direction of the plane
	Cart3D three_d;

	// dimensions needed to create a two-dimensional projection
	// of the top surface
	Cart2D two_d;

	// map from a two-dimensional Cartesian coordinate to the final
	// id of an active element in the grid, or -1 if nothing is assigned
	// this vector is first valid after create_elements() have been done
	vector <int> elms;

	// map from a two-dimensional Cartesian node coordinate to the final
	// id of active nodes in the grid, or -1 if nothing is assigned.
	// this vector is first valid after create_nodes() have been done
	vector <int> nodes;

	TopSurfBuilder (const UnstructuredGrid& from, TopSurf& into)
		// link to the fine grid for the duration of the construction
		: fine_grid (from)

		// allocate memory for the grid. it is initially empty
		, ts (into)

		// extract dimensions from the source grid
		, three_d (fine_grid)
		, two_d (three_d.project ()) {

		// check that the fine grid contains structured information;
		// this is essential to mapping cells to columns
		if (!fine_grid.global_cell) {
			throw OPM_EXC ("Find grid is not (logically) structured");
		}

		// create frame of the new top surface
		create_dimensions ();

		// identify active columns in the grid
		create_elements ();

		// identify active points in the grid
		create_nodes ();
	}

	// various stages of the build process, supposed to be called in
	// this order. (I have separated them into separate procedures to
	// make it more obvious what parts that needs to be shared between
	// them)
private:
	void create_dimensions () {
		// we are going to create two-dimensional grid
		ts.dimensions = 2;
		ts.cartdims[0] = two_d.ni;
		ts.cartdims[1] = two_d.nj;
		ts.cartdims[2] = 1;
	}

	void create_elements() {
		// statistics of the deepest and highest active k-index of
		// each column in the grid. to know each index into the column,
		// we only need to know the deepest k and the count; the highest
		// is sampled to do consistency checks afterwards
		const int num_cols = two_d.num_elems ();

		// assume initially that there are no active elements in each column
		vector <int> act_cnt (num_cols, 0);
		ts.number_of_cells = 0;

		// initialize these to values that are surely out of range, so that
		// the first invocation of min or max always set the value. we use
		// this to detect whether anything was written later on. since the
		// numbering of the grid starts at the bottom, then the deepest cell
		// has the *lowest* k-index, thus we need a value greater than all
		vector <int> deep_k (num_cols, INT_MAX);
		vector <int> high_k (num_cols, INT_MIN);

		// loop once through the fine grid to gather statistics of the
		// size of the surface so we know what to allocate
		for (int fine_elem = 0; fine_elem != fine_grid.number_of_cells; ++fine_elem) {
			// get the cartesian index for this cell; this is the cell
			// number in a grid that also includes the inactive cells
			const Cart3D::elem_t cart_ndx = fine_grid.global_cell [fine_elem];

			// deconstruct the cartesian index into (i,j,k) constituents;
			// the i-index moves fastest, as this is Fortran-indexing
			const Coord3D ijk = three_d.coord (cart_ndx);

			// figure out which column this item belongs to (in 2D)
			const Cart2D::elem_t col = two_d.cart_ndx (ijk);

			// update the statistics for this column; 'deepest' is the lowest
			// k-index seen so far, 'highest' is the highest (ehm)
			deep_k[col] = min (deep_k[col], ijk.k);
			high_k[col] = max (high_k[col], ijk.k);

			// we have seen an element in this column; it becomes active. only
			// columns with active cells will get active elements in the surface
			// grid. if this cell wasn't marked as active before we can check it
			// off now
			if (!act_cnt[col]++) {
				ts.number_of_cells++;
			}
		}

		// check that we have a continuous range of elements in each column;
		// this must be the case to assume that the entire column can be merged
		for (int col = 0; col < num_cols; ++col) {
			if (act_cnt[col]) {
				if (deep_k[col] + act_cnt[col] - 1 != high_k[col]) {
					const Coord2D coord = two_d.coord (col);
					throw OPM_EXC ("Non-continuous column at (%d, %d)", coord.i, coord.j);
				}
			}
		}

		// allocate memory needed to hold the elements in the grid structure
		// if we throw an exception at some point, the destructor of the TopSurf
		// will take care of deallocating this memory for us
		ts.global_cell = new int [ts.number_of_cells];
		ts.col_cellpos = new int [ts.number_of_cells+1];

		// there is no elements before the first column, so this number is
		// always zero. it is written to the array to maintain a straight code
		// path in calculations.
		ts.col_cellpos[0] = 0;

		// now we know enough to start assigning ids to active columns.if
		// memory is a real shortage, we could reuse the act_cnt array for this.
		elms.resize (num_cols, -1);

		// loop through the grid and assign an id for all columns that have
		// active elements
		int elem_id = 0;
		for (int col = 0; col < num_cols; ++col) {
			if (act_cnt[col]) {
				elms[col] = elem_id;

				// dual pointer that maps the other way; what is the structured
				// index of this element (flattened into an integer)
				ts.global_cell[elem_id] = col;

				// note the number of elements there now are before the next column;
				// in addition to all the previous ones, our elements are now added
				ts.col_cellpos[elem_id+1] = ts.col_cellpos[elem_id] + act_cnt[col];
				elem_id++;
			}
		}

		// now write indices from the fine grid into the column map of the surface
		// we end up with a list of element that are in each column
		ts.col_cells = new int [fine_grid.number_of_cells];
		for (int cell = 0; cell < fine_grid.number_of_cells; ++cell) {
			// get the Cartesian index for this element
			const Cart3D::elem_t cart_ndx = fine_grid.global_cell[cell];
			const Coord3D ijk = three_d.coord (cart_ndx);

			// get the id of the column in which this element now belongs
			const Cart2D::elem_t col = two_d.cart_ndx (ijk);
			const int elem_id = elms[col];

			// start of the list of elements for this particular column
			const int segment = ts.col_cellpos[elem_id];

			// since there is supposed to be a continuous range of elements in
			// each column, we can calculate the relative position in the list
			// based on the k part of the coordinate.
			const int offset = ijk.k - deep_k[col];

			// write the fine grid cell number in the column list; since we
			// have calculated the position based on depth, the list will be
			// sorted downwards up when we are done
			ts.col_cells[segment + offset] = cell;
		}

		// these members will be filled by computeGeometry, but needs valid
		// memory to work with
		ts.cell_volumes = new double [ts.number_of_cells];
		ts.cell_centroids = new double [ts.dimensions * ts.number_of_cells];
	}

	void create_nodes () {
		// construct a dual Cartesian grid consisting of the points
		const int num_nodes = two_d.num_nodes ();

		// vectors which will hold the coordinates for each active point.
		// at first we sum all the points, then we divide by the count to
		// get the average. as long as the count is zero, there is no
		// registered active point at this location
		vector <double> x (num_nodes, 0.);
		vector <double> y (num_nodes, 0.);
		vector <int> cnt (num_nodes, 0);

		// number of nodes needed in the top surface
		int active_nodes = 0;

		// tag of the top side in a cell; we're looking for this
		const int top_tag = Side3D (Dim3D::Z, Dir::INC).facetag ();

		// initial corner value. this could really be anything, since
		// we expect all the fields to be overwritten.
		const Corn3D blank (Dir::DEC, Dir::DEC, Dir::DEC);

		// this map holds the classification of each node locally for the
		// element being currently processed. we reuse the map to avoid
		// needless memory allocations.
		typedef map <int, Corn3D> cls_t;
		cls_t classifier;

		// loop through all active cells in the top surface
		for (int col = 0; col < ts.number_of_cells; ++col) {
			// get the highest element in this column; since we have them
			// sorted by k-index this should be the last item in the
			// extended column info, or one before the next column starts
			const Cart3D::elem_t top_cell =
					ts.col_cells [ts.col_cellpos[col+1] - 1];

			// start afresh whenever we start working on a new element
			classifier.clear ();
			int top_face = -1;

			// loop through all the faces of the top element
			for (int face_pos = fine_grid.cell_facepos[top_cell];
					 face_pos != fine_grid.cell_facepos[top_cell+1];
					 ++face_pos) {

				// get the (normal) dimension and direction of this face
				const int this_tag = fine_grid.cell_facetag[face_pos];
				Side3D s = Side3D::from_tag (this_tag);

				// identifier of the face, which is the index in the next arary
				const int face = fine_grid.cell_faces[face_pos];

				// remember it if we've found the top face
				if (this_tag == top_tag) {
					if (top_face != -1) {
						throw OPM_EXC ("More than one top face in element %d", top_cell);
					}
					top_face = face;
				}

				// loop through all nodes in this face, adding them to the
				// classifier. when we are through with all the faces, we have
				// found in which corner a node is, defined by a direction in
				// each of the three dimensions
				for (int node_pos = fine_grid.face_nodepos[face];
						 node_pos != fine_grid.face_nodepos[face+1];
						 ++node_pos) {
					const int node = fine_grid.face_nodes[node_pos];

					// locate pointer to data record ("iterator" in stl parlance)
					// for this node, if it is already there. otherwise, just start
					// out with some blank data (which eventually will get overwritten)
					cls_t::iterator ptr = classifier.find (node);
					Corn3D prev (ptr == classifier.end () ? blank : ptr->second);

					// update the dimension in which this face is pointing
					if (ptr != classifier.end ()) {
						classifier.erase (ptr);
					}
					classifier.insert (ptr, make_pair (node, prev.pivot (s.dim, s.dir)));
				}

				// after this loop, we have a map of each node local to the element,
				// classified into in which corner it is located (it cannot be in
				// both directions in the same dimension -- then it would have to
				// belong to two opposite faces, unless the grid is degenerated)
			}

			// cannot handle degenerate grids without top face properly
			if (top_face == -1) {
				throw OPM_EXC ("No top face in cell %d", top_cell);
			}

			// get the Cartesian ij coordinate of this cell
			const Coord2D ij = two_d.coord (ts.global_cell [top_cell]);

			// loop through all the nodes of the top face, and write their position
			// into the corresponding two-d node. this has to be done separately
			// after we have classified *all* the nodes of the element, in order for
			// the corner values to be set correctly, i.e. we cannot merge this into
			// the loop above.
			for (int node_pos = fine_grid.face_nodepos[top_face];
					 node_pos != fine_grid.face_nodepos[top_face+1];
					 ++node_pos) {
				const int node = fine_grid.face_nodes[node_pos];

				// get which corner this node has; this returns a three-dimensional
				// corner, but by using the base class part of it we automatically
				// project it to a flat surface
				cls_t::iterator ptr = classifier.find (node);
				const Corn3D corn (ptr->second);

				// get the structured index for this particular corner
				const Cart2D::node_t cart_node = two_d.node_ndx(ij, corn);

				// add these coordinates to the average position for this junction
				// if we activate a corner, then add it to the total count
				x[cart_node] += fine_grid.node_coordinates[node+0];
				y[cart_node] += fine_grid.node_coordinates[node+1];
				if (!cnt[cart_node]++) {
					++active_nodes;
				}
			}
		}

		// after this loop we the accumulated coordinates for each of the
		// corners that are part of active elements (the nodes that are
		// needed in the top surface)

		// assign identifiers and find average coordinate for each point
		ts.number_of_nodes = active_nodes;
		ts.node_coordinates = new double [active_nodes * Dim2D::COUNT];
		nodes.resize (num_nodes, -1);
		int next_node_id = 0;
		for (int cart_node = 0; cart_node != num_nodes; ++cart_node) {
			if (cnt[cart_node]) {
				nodes[cart_node] = next_node_id;
				const int start = Dim2D::COUNT * next_node_id;
				ts.node_coordinates[start+0] = x[cart_node] / cnt[cart_node];
				ts.node_coordinates[start+1] = y[cart_node] / cnt[cart_node];
				++next_node_id;
			}
		}
	}
};

TopSurf*
TopSurf::create (const UnstructuredGrid& fine_grid) {
	// I *know* that we are not supposed to use auto_ptr anymore, but
	// it works in both C++98 and C++11 and is the only common way to
	// transfer the object out of the smart pointer afterwards
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
	auto_ptr <TopSurf> ts (new TopSurf);
#pragma GCC diagnostic pop

	// outsource the entire construction to a builder object
	TopSurfBuilder (fine_grid, *(ts.get ()));

	// client owns pointer to constructed grid from this point
	return ts.release ();
}

TopSurf::TopSurf ()
	: col_cells (0)
	, col_cellpos (0) {
	// zero initialize all members that come from UnstructuredGrid
	// since that struct is a C struct, it doesn't have a ctor
	dimensions = 0;
	number_of_cells = 0;
	number_of_faces = 0;
	number_of_nodes = 0;
	face_nodes = 0;
	face_nodepos = 0;
	face_cells = 0;
	cell_faces = 0;
	cell_facepos = 0;
	node_coordinates = 0;
	face_centroids = 0;
	face_areas = 0;
	face_normals = 0;
	cell_centroids = 0;
	cell_volumes = 0;
	global_cell = 0;
	cartdims[0] = 0;
	cartdims[1] = 0;
	cartdims[2] = 0;
	cell_facetag = 0;
}

TopSurf::~TopSurf () {
	// deallocate memory that may have been created. if the dtor is
	// called from throwing an exception, the members should be zero
	// initialized, so it's OK to send them to delete.
	delete [] face_nodes;
	delete [] face_nodepos;
	delete [] face_cells;
	delete [] cell_faces;
	delete [] cell_facepos;
	delete [] node_coordinates;
	delete [] face_centroids;
	delete [] face_areas;
	delete [] face_normals;
	delete [] cell_volumes;
	delete [] global_cell;
	delete [] cell_facetag;
	// these are the extra members that are TopSurf specific
	delete [] col_cells;
	delete [] col_cellpos;
}
