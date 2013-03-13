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
