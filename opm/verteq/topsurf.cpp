// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#include <opm/verteq/topsurf.hpp>
#include <opm/core/utility/exc.hpp>
#include <boost/range/iterator_range.hpp>
#include <memory> // auto_ptr

using namespace boost;
using namespace Opm;
using namespace std;

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

	TopSurfBuilder (const UnstructuredGrid& from, TopSurf& into)
		// link to the fine grid for the duration of the construction
		: fine_grid (from)

		// allocate memory for the grid. it is initially empty
		, ts (into) {

		// check that the fine grid contains structured information;
		// this is essential to mapping cells to columns
		if (!fine_grid.global_cell) {
			throw OPM_EXC ("Find grid is not (logically) structured");
		}
	}

	// various stages of the build process, supposed to be called in
	// this order. (I have separated them into separate procedures to
	// make it more obvious what parts that needs to be shared between
	// them)
private:
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

iterator_range <int*>
TopSurf::column (int ndx_2d) {
	// skip all the cells that belongs to columns before us,
	// then we get to the start of our own column
	int* begin_addr = &this->col_cells [col_cellpos [ndx_2d]];

	// stop when we arrive at the start of the next column
	int* end_addr = &this->col_cells [col_cellpos [ndx_2d + 1]];

	// return an iterator over this
	return iterator_range <int*> (begin_addr, end_addr);
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
