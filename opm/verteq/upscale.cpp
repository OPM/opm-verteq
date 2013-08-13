#include <opm/verteq/upscale.hpp>
#include <opm/verteq/utility/runlen.hpp>
#include <cmath> // floor

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

double
VertEqUpscaler::sum (
		int col,
		const double *val) const {

	// index into the fine grid for all columns
	const rlw_int col_cells (ts.number_of_cells, ts.col_cellpos, ts.col_cells);

	// get the indices for this particular column
	const int* fine_ndx = col_cells[col];

	// loop through each block in the column, accumulating as we go
	double sum = 0.;
	for (int row = 0; row < col_cells.size (col); ++row) {
		sum += val[fine_ndx[row]];
	}

	return sum;
}

void
VertEqUpscaler::wgt_dpt (
		int col,
		const double* val,
		double* res) const {

	// get the weights for this particular column
	const rlw_double dz = rlw_double (ts.number_of_cells, ts.col_cellpos, ts.dz);
	const double* dz_col = dz[col];

	// running total
	double accum = 0.;

	// divisor for this particular column. this should be larger than zero
	// since the column otherwise wouldn't be active
	double H = ts.h_tot[col];

	// loop through each of the blocks from the top
	for (int row = 0; row < dz.size(col); ++row) {
		// this is \int_{h}^{\zeta_T} f
		accum += val[row] * dz_col[row];

		// write the average to the resulting array
		res[row] = accum / H;
	}
}

double
VertEqUpscaler::dpt_avg (
		int col,
		const double* val) const {

	// get the weights for this particular column
	const rlw_double dz = rlw_double (ts.number_of_cells, ts.col_cellpos, ts.dz);
	const double* dz_col = dz[col];

	// running total
	double accum = 0.;

	// multiply each row with its depth
	for (int row = 0; row < dz.size(col); ++row) {
		accum += val[row] * dz_col[row];
	}

	// we already have stored the total height of the column; no need to
	// keep a counter for that too
	const double avg = accum / ts.h_tot[col];
	return avg;
}

int
VertEqUpscaler::num_rows (
    int col) const {

	// use this helper object to query about the size of the column
	// (the compiler should be able to optimize most of it away)
	const rlw_int pos (ts.number_of_cells, ts.col_cellpos, ts.col_cells);
	return pos.size (col);
}

Elevation
VertEqUpscaler::bottom (
    int col) const {

	// simply initialize to skip *all* blocks in that column
	return Elevation (num_rows (col), 0.);
}

double
VertEqUpscaler::eval (
		int col,
		const double* dpt,
    const Elevation zeta) const {

	// number of whole blocks to include
	const int row = zeta.block ();

	// if any blocks before, take those values. unfortunately the implied
	// zero value in front of the array cause a branch; hopefully branch
	// prediction in the CPU will be able to alliviate some of that cost.
	// if we add an explicit zero, then the memory sizes to store the values
	// are not the same as in the top grid, and we'll have to adjust indices
	// all the time.
	const double before = row == 0 ? 0. : dpt[row-1];

	// then we add the fractional value of just this block. it is just the
	// last block that we want the fraction for
	return before + (dpt[row] - before) * zeta.fraction ();
}

Elevation
VertEqUpscaler::find (
		int col,
		const double* dpt,
		const double target) const {

	// use interpolation search to find the proper block for this
	// (relative) depth. under the assumption that all individual
	// heights are of approx. the same order, this should have time
	// complexity O(log log N), with good locality. however, the
	// constant factor is higher than a binary search due to the
	// floating point arithmetic in computing the guess.

	// we start searching from the first value; the top_hgt is the
	// height at the *start* of the block (unlike what is stored in
	// the up_bnd array)
	int top_ndx = 0;
	double top_val = 0.; // top_hgt <= target

	// bottom of the searching scope. bot_hgt is the height of the
	// *end* of the block, like what is in up_bnd
	int bot_ndx = num_rows (col) - 1;
	double bot_val = dpt[bot_ndx]; // target <=bot_val

	// search until we have a solution; as the search space get tighter
	// and tighter, we must eventually end up with a solution if the
	// saturation was between 0. and 1 initially.
	for (;;) {

		// based on the fraction of the interval the target height is,
		// guess at the index assuming that every block has equal height
		const double frac = (target - top_val) / (bot_val - top_val);
		const int cur_ndx = top_ndx + static_cast <int> (
		                    std::floor ((bot_ndx - top_ndx + 1) * frac));

		// get the brackets of this block; the weigted depth is the upper
		// bound of the integral for each block. unfortunately we don't have
		// lower bounds stored in an array, so we must have a conditional
		const double cur_bot = dpt[cur_ndx];
		const double cur_top = cur_ndx == 0 ? 0. : dpt[cur_ndx - 1];

		// divide the search space into three: the current block, everything
		// before and everything after. if we are in the current bracket,
		// then return the result, otherwise search further in one of the
		// two subspaces.
		if (cur_bot < target) {
			top_ndx = cur_ndx + 1;
			top_val = cur_bot; // top_val still < target
		}
		else if (target < cur_top) {
			bot_ndx = cur_ndx - 1;
			bot_val = cur_top; // bot_val still > target
		}
		else {
			// get the fraction of this block
			const double cur_frac = (target - cur_top) / (cur_bot - cur_top);
			return Elevation (cur_ndx, cur_frac);
		}
	}
}
