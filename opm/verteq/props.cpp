#include <opm/verteq/props.hpp>
#include <opm/verteq/topsurf.hpp>
#include <opm/verteq/upscale.hpp>
#include <opm/verteq/utility/exc.hpp>
#include <algorithm> // fill
#include <memory> // auto_ptr
using namespace Opm;
using namespace std;

struct VertEqPropsImpl : public VertEqProps {
	/// Get the underlaying fluid information from here
	const IncompPropertiesInterface& fp;

	/// Get the grid information from here
	const TopSurf& ts;

	// constants to avoid a bunch of "magic" values in the code
	static const int TWO_DIMS   = 2;
	static const int THREE_DIMS = 3;

	// size of the permeability matrices (in numbers)
	static const int PERM_MATRIX_2D = TWO_DIMS * TWO_DIMS;
	static const int PERM_MATRIX_3D = THREE_DIMS * THREE_DIMS;

	// offsets when indexing into the permeability matrix
	static const int KXX_OFS_3D = 0 * THREE_DIMS + 0; // (x, x), x = 0
	static const int KXY_OFS_3D = 0 * THREE_DIMS + 1; // (x, y), x = 0, y = 1
	static const int KYY_OFS_3D = 1 * THREE_DIMS + 1; // (y, y), y = 1

	static const int KXX_OFS_2D = 0 * TWO_DIMS + 0; // (x, x), x = 0
	static const int KXY_OFS_2D = 0 * TWO_DIMS + 1; // (x, y), x = 0, y = 1
	static const int KYX_OFS_2D = 1 * TWO_DIMS + 0; // (y, x), x = 0, y = 1
	static const int KYY_OFS_2D = 1 * TWO_DIMS + 1; // (y, y), y = 1

	/// Helper object to do averaging
	const VertEqUpscaler up;

	/// Upscaled porosity; this is \Phi in the papers
	vector <double> upscaled_poro;

	/// Upscaled permeability; this is K in the papers
	vector <double> upscaled_absperm;

	VertEqPropsImpl (const IncompPropertiesInterface& fineProps,
	                 const TopSurf& topSurf)
		: fp (fineProps)
		, ts (topSurf)
		, up (ts) {

		// allocate memory to store results for faster lookup later
		upscaled_poro.resize (ts.number_of_cells);

		// buffers that holds intermediate values for each column;
		// pre-allocate to avoid doing that inside the loop
		vector <double> poro (ts.max_vert_res, 0.); // porosity
		vector <double> kxx (ts.max_vert_res, 0.);  // abs.perm.
		vector <double> kxy (ts.max_vert_res, 0.);
		vector <double> kyy (ts.max_vert_res, 0.);

		// pointer to all porosities in the fine grid
		const double* fine_poro = fp.porosity ();
		const double* fine_perm = fp.permeability ();

		// upscale each column separately
		for (int col = 0; col < ts.number_of_cells; ++col) {
			// retrieve the fine porosities for this column only
			up.gather (col, &poro[0], fine_poro, 1, 0);

			// compute the depth-averaged value and store
			upscaled_poro[col] = up.dpt_avg (col, &poro[0]);

			// retrieve the fine abs. perm. for this column only
			up.gather (col, &kxx[0], fine_perm, PERM_MATRIX_3D, KXX_OFS_3D);
			up.gather (col, &kxy[0], fine_perm, PERM_MATRIX_3D, KXY_OFS_3D);
			up.gather (col, &kyy[0], fine_perm, PERM_MATRIX_3D, KYY_OFS_3D);

			// compute upscaled values for each dimension separately
			const double up_kxx = up.dpt_avg (col, &kxx[0]);
			const double up_kxy = up.dpt_avg (col, &kxy[0]);
			const double up_kyy = up.dpt_avg (col, &kyy[0]);

			// store back into the interleaved format required by the 2D
			// simulator code (fetching a tensor at the time, probably)
			// notice that we take advantage of the tensor being symmetric
			// at the third line below
			upscaled_absperm[PERM_MATRIX_2D * col + KXX_OFS_2D] = up_kxx;
			upscaled_absperm[PERM_MATRIX_2D * col + KXY_OFS_2D] = up_kxy;
			upscaled_absperm[PERM_MATRIX_2D * col + KYX_OFS_2D] = up_kxy;
			upscaled_absperm[PERM_MATRIX_2D * col + KYY_OFS_2D] = up_kyy;
		}
	}

	/* rock properties; use volume-weighted averages */
	virtual int numDimensions () const {
		// the upscaled grid is always dimensionally reduced
		return TWO_DIMS;
	}

	virtual int numCells () const {
		// we'll provide on value for each column in the upscaled grid
		return ts.number_of_cells;
	}

	virtual const double* porosity () const {
		// calculated in the constructor; since we must return a full
		// array there isn't anything to save by calculating on the fly
		// (accessing the data like this is supposed to be safe as long
		// as the container "lives")
		return &upscaled_poro[0];
	}

	virtual const double* permeability () const {
		return &upscaled_absperm[0];
	}

	/* fluid properties; these don't change when upscaling */
	virtual int numPhases () const {
		return fp.numPhases ();
	}

	virtual const double* viscosity () const {
		return fp.viscosity ();
	}

	virtual const double* density () const {
		return fp.density ();
	}

	virtual const double* surfaceDensity () const {
		return fp.surfaceDensity ();
	}

	/* hydrological (unsaturated zone) properties */
	virtual void relperm (const int n,
	                      const double *s,
	                      const int *cells,
	                      double *kr,
	                      double *dkrds) const {
		throw OPM_EXC ("Not implemented yet");
	}

	virtual void capPress (const int n,
	                       const double *s,
	                       const int *cells,
	                       double *pc,
	                       double *dpcds) const {
		throw OPM_EXC ("Not implemented yet");
	}

	virtual void satRange (const int n,
	                       const int *cells,
	                       double *smin,
	                       double *smax) const {
		// saturation is just another name for "how much of the column
		// is filled", so every range from nothing to completely filled
		// are valid. even though there is residual water/gas in each
		// block, this is not seen from the 2D code
		const int np = n * numPhases ();
		fill (smin, smin + np, 0.);
		fill (smax, smax + np, 1.);
	}
};

VertEqProps*
VertEqProps::create (const IncompPropertiesInterface& fineProps,
                     const TopSurf& topSurf) {
	// construct real object which contains all the implementation details
	auto_ptr <VertEqProps> props (new VertEqPropsImpl (fineProps, topSurf));

	// client owns pointer to constructed fluid object from this point
	return props.release ();
}
