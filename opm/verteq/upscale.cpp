#include <opm/verteq/upscale.hpp>
#include <opm/verteq/utility/runlen.hpp>

using namespace Opm;

void
VertEqUpscaler::gather (
		int col,
		double* buf,
		const double* data,
		int stride,
		int offset) const {

	// index into the fine grid for all columns
	const rlw_int col_cells (ts.number_of_cells, ts.col_cellpos, ts.col_cells);

	// get the indices for this particular column
	const int* fine_ndx = col_cells[col];

	// loop through each block in the column and fetch the property
	for (int row = 0; row < col_cells.size (col); ++row) {
		// index in the fine grid for this block
		const int block_ndx = fine_ndx[row];

		// calculate position in the data array
		const int pos = block_ndx * stride + offset;

		// move the data
		buf[row] = data[pos];
	}
}
