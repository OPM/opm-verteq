// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#include <opm/verteq/topsurf.hpp>
#include <boost/range/iterator_range.hpp>
#include <memory> // auto_ptr

using namespace boost;
using namespace Opm;
using namespace std;

TopSurf*
TopSurf::create (const UnstructuredGrid& fine) {
	// allocate memory for the grid. it is initially empty
	auto_ptr <TopSurf> ts = auto_ptr <TopSurf> (new TopSurf);

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
