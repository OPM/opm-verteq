#include <opm/verteq/props.hpp>
#include <opm/verteq/topsurf.hpp>
#include <opm/verteq/upscale.hpp>
#include <opm/verteq/utility/exc.hpp>
#include <opm/verteq/utility/runlen.hpp>
#include <algorithm> // fill
#include <cmath> // sqrt
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

	// we assume this ordering of the phases in arrays
	static const int GAS = 0;
	static const int WAT = 1;

	/// Helper object to do averaging
	const VertEqUpscaler up;

	/// Upscaled porosity; this is \Phi in the papers
	vector <double> upscaled_poro;

	/// Upscaled permeability; this is K in the papers
	vector <double> upscaled_absperm;

	/// Volume fractions of gas phase, used in averaging
	RunLenData <double> res_gas_vol; // \phi S_{n,r}
	RunLenData <double> mob_mix_vol; // \phi (1 - S_{w,r} - S_{n,r})
	RunLenData <double> res_wat_vol; // \phi (1 - S_{w,r})

	/// Volume-of-gas-phase-fraction-weighted depths-fractions
	RunLenData <double> res_gas_dpt; // int_{h}^{\zeta_T} \phi S_{n,r} dz
	RunLenData <double> mob_mix_dpt; // int_{h}^{\zeta_T} \phi (1 - S_{w,r} - S_{n,r} dz
	RunLenData <double> res_wat_dpt; // int_{h}^{\zeta_T} \phi (1 - S_{w,r}) dz

	// we need to keep track of where the plume has been and deposited
	// residual CO2. however, finding the interface is non-trivial and
	// should only be done if we actually see a new maximum of the
	// saturation. this array contains the trigger point for recalc.
	vector <double> max_gas_sat;      // S_{g,max}
	vector <Elevation> max_gas_elev;  // \zeta_R

	virtual void upd_res_sat (const double* snap) {
		// cache this here outside of the loop
		const int num_phases = numPhases ();

		// update saturation for each column
		for (int col = 0; col < ts.number_of_cells; ++col) {
			// current CO2 saturation
			const double cur_sat = snap[col * num_phases + GAS];

			// has it increased? is there more of the plume in this column?
			check_res_sat (col, cur_sat);
		}
	}

	void check_res_sat (int col, double cur_sat) {
		if (cur_sat > max_gas_sat[col]) {
			// recalculate discretized elevation
			max_gas_elev[col] = res_elev (col, cur_sat);

			// update stored saturation so we test correctly next time
			max_gas_sat[col] = cur_sat;
		}
	}

	/**
	 * Find the elevation of the residual CO2 in this column based on the
	 * maximum upscaled CO2 saturation.
	 *
	 * This is done by solving this equation for \zeta_R:
	 *
	 * H \Phi S_{g,max} = \int_{\zeta_R}^{\zeta_T} \phi (1 - s_{w,r}) dz
	 *
	 * using precalculated values for the integral.
	 */
	Elevation res_elev(const int col, const double max_sat) {
		// right-hand side of the equation (apart from H, which is divided
		// in the averaging operator stored)
		const double max_vol = upscaled_poro[col] * max_sat;

		// find the elevation which makes the integral have this value
		const Elevation zeta_r = up.find (col, res_wat_dpt[col], max_vol);
		return zeta_r;
	}

	/**
	 * Find the elevation (as a table index) for the current interface
	 * between the plume of CO2 and the brine layer below, given an
	 * upscaled saturation of CO2 (which is correlated to the height)
	 *
	 * This is done by solving this equation for \zetaM:
	 *
	 * H \Phi S_g =    \int_{\zeta_R}^{\zeta_T} \phi s_{g,r} dz
	 *               + \int_{\zeta_M}^{\zeta_T} \phi (1-s_{w,r}-s_{g,r} dz
	 *
	 * using precalculated values for the integrals, and a stored table
	 * index for the extent of the
	 *
	 * This should be done *after* the maximum saturation is updated,
	 * so we have \zeta_R^(t) and not \zeta_R^(t-1) (which may be higher
	 * up than the current interface!)
	 */
	Elevation intf_elev (const int col, const double gas_sat) const {
		// check to make sure that the simulator updates the correct values
		// this duplicates the efforts of the callback, but this is only called
		// whenever we need the rel.perm.
		//check_res_sat (col, gas_sat);

		// the first term is \Phi * S_g representing the volume of CO2, the
		// second is the integral int_{\zeta_R}^{\zeta_T} \phi s_{g,r} dz,
		// representing the volume of residual CO2; the remainder becomes
		// the mobile CO2 volume
		const double gas_vol = upscaled_poro[col] * gas_sat // \Phi * S_g
		                     + up.eval (col, res_gas_dpt[col], max_gas_elev[col]);

		// lookup to find the height that gives this mobile volume
		const Elevation zeta_M = up.find (col, mob_mix_dpt[col], gas_vol);
		return zeta_M;
	}

	/**
	 * Magnitude of (lateral) permeability tensor used in weighting.
	 *
	 * Since the upscaled relative permeabilities is a weighting of
	 * the various fine-scale absolute ones, they will end up as tensors
	 * as well. However, our simulator code needs a scalar to measure
	 * the "weight" of an individual tensor.
	 *
	 * We choose the Frobenius norm as the contraction operation on the
	 * tensor.
	 */
	double magnitude (double kxx, double kxy, double kyy) {
		return sqrt (kxx*kxx + 2*kxy*kxy + kyy*kyy);
	}

	// weighted rel.perm. for CO2 when residual brine is present, and the
	// depth for each block further weighted with this.
	RunLenData <double> prm_gas;      // K^{-1} k_|| k_{g,r} (1-s_{w,r})
	RunLenData <double> prm_gas_int;  // 1/H \int_h^{\zeta_T} above dz

	// weighted rel.perm. for residual part of brine
	RunLenData <double> prm_res;      // K^{-1} k_|| 1 - k_{w,r} (s_{g,r})
	RunLenData <double> prm_res_int;  // 1/H \int_h^{\zeta_T} above dz

	// weighted rel.perm. for brine when residual CO2 is present
	RunLenData <double> prm_wat;      // K^{-1} k_|| k_{w,r} (s_{g,r})
	RunLenData <double> prm_wat_int;  // 1/H \int_h^{\zeta_T} above dz

	// gravity in the z-direction; \nabla z \cdot \mathbf{g}
	const double gravity;

	VertEqPropsImpl (const IncompPropertiesInterface& fineProps,
	                 const TopSurf& topSurf,
	                 const double* grav_vec)
		: fp (fineProps)
		, ts (topSurf)
		, up (ts)
		, res_gas_vol (ts.number_of_cells, ts.col_cellpos)
		, mob_mix_vol (ts.number_of_cells, ts.col_cellpos)
		, res_wat_vol (ts.number_of_cells, ts.col_cellpos)
		, res_gas_dpt (ts.number_of_cells, ts.col_cellpos)
		, mob_mix_dpt (ts.number_of_cells, ts.col_cellpos)
		, res_wat_dpt (ts.number_of_cells, ts.col_cellpos)

		// assume that there is no initial plume; first notification will
		// trigger an update of all columns where there actually is CO2
		, max_gas_sat (ts.number_of_cells, 0.)

		// this is the elevation that corresponds to no CO2 sat.
		, max_gas_elev (ts.number_of_cells, Elevation (0, 0.))

		, prm_gas (ts.number_of_cells, ts.col_cellpos)
		, prm_gas_int (ts.number_of_cells, ts.col_cellpos)
		, prm_res (ts.number_of_cells, ts.col_cellpos)
		, prm_res_int (ts.number_of_cells, ts.col_cellpos)
		, prm_wat (ts.number_of_cells, ts.col_cellpos)
		, prm_wat_int (ts.number_of_cells, ts.col_cellpos)
		, gravity (grav_vec[THREE_DIMS - 1])	{

		// allocate memory to store results for faster lookup later
		upscaled_poro.resize (ts.number_of_cells);

		// number of phases (should be 2)
		const int num_phases = fp.numPhases ();

		// buffers that holds intermediate values for each column;
		// pre-allocate to avoid doing that inside the loop
		vector <double> poro (ts.max_vert_res, 0.); // porosity
		vector <double> kxx (ts.max_vert_res, 0.);  // abs.perm.
		vector <double> kxy (ts.max_vert_res, 0.);
		vector <double> kyy (ts.max_vert_res, 0.);
		vector <double> sgr   (ts.max_vert_res * num_phases, 0.); // residual CO2
		vector <double> l_swr (ts.max_vert_res * num_phases, 0.); // 1 - residual brine
		vector <double> lkl (ts.max_vert_res); // magnitude of abs.perm.; k_||

		// saturations and rel.perms. of each phase, assuming maximum filling of...
		vector <double> wat_sat (ts.max_vert_res * num_phases, 0.); // brine; res. CO2
		vector <double> gas_sat (ts.max_vert_res * num_phases, 0.); // CO2; res. brine
		vector <double> wat_mob (ts.max_vert_res * num_phases, 0.); // k_r(S_c=S_{c,r})
		vector <double> gas_mob (ts.max_vert_res * num_phases, 0.); // k_r(S_c=1-S_{b,r})

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

			// contract each fine perm. to a scalar, used for weight later
			for (int row = 0; row < up.num_rows (col); ++row) {
				lkl[row] = magnitude (kxx[row], kxy[row], kyy[row]);
			}

			// we only need the relative weight, so get the depth-averaged
			// total weight, which we'll use to scale the weights below
			const double tot_lkl = up.dpt_avg (col, &lkl[0]); // 1/K^{-1}

			// query the fine properties for the residual saturations;
			// notice that we implicitly get the brine saturation as the maximum
			// allowable co2 saturation; now we've got the values we need, but
			// only every other item (due to that both phases are stored)
			const rlw_int col_cells (ts.number_of_cells, ts.col_cellpos, ts.col_cells);
			fp.satRange (col_cells.size (col), col_cells[col], &sgr[0], &l_swr[0]);

			// cache pointers to this particular column to avoid recomputing
			// the starting point for each and every item
			double* res_gas_col = res_gas_vol[col];
			double* mob_mix_col = mob_mix_vol[col];
			double* res_wat_col = res_wat_vol[col];

			for (int row = 0; row < col_cells.size (col); ++row) {
				// multiply with num_phases because the saturations for *both*
				// phases are store consequtively (as a record); we only need
				// the residuals framed as co2 saturations
				const double sgr_ = sgr[row * num_phases + GAS];
				const double l_swr_ = l_swr[row * num_phases + GAS];

				// portions of the block that are filled with: residual co2,
				// mobile fluid and residual brine, respectively
				res_gas_col[row] = poro[row] * sgr_;            // \phi*S_{n,r}
				mob_mix_col[row] = poro[row] * (l_swr_ - sgr_); // \phi*(1-S_{w,r}-S_{n_r})
				res_wat_col[row] = poro[row] * l_swr_;          // \phi*(1-S_{w,r}
			}

			// weight the relative depth factor (how close are we towards a
			// completely filled column) with the volume portions
			up.wgt_dpt (col, &res_gas_col[0], &res_gas_dpt[col][0]);
			up.wgt_dpt (col, &mob_mix_col[0], &mob_mix_dpt[col][0]);
			up.wgt_dpt (col, &res_gas_col[0], &res_wat_dpt[col][0]);

			// now, when we queried the saturation ranges, we got back the min.
			// and max. sat., and when there is min. of one, then there should
			// be max. of the other; however, these data are in different arrays!
			// cross-pick such that we get (min CO2, max brine), (max CO2, min brine)
			// instead of (min CO2, min brine), (max CO2, max brine). this code
			// has no other effect than to satisfy the ordering of items required
			// for the relperm() call
			for (int row = 0; row < col_cells.size (col); ++row) {
				wat_sat[row * num_phases + GAS] = sgr[row * num_phases + GAS];
				wat_sat[row * num_phases + WAT] = l_swr[row * num_phases + WAT];
				gas_sat[row * num_phases + GAS] = l_swr[row * num_phases + GAS];
				gas_sat[row * num_phases + WAT] = sgr[row * num_phases + WAT];
			}

			// get rel.perm. for those cases where one phase is (maximally) mobile
			// and the other one is immobile (at residual saturation); we get back
			// rel.perm. for both phases, although only one of them is of interest
			// for us (the other one should be zero). we have no interest in the
			// derivative of the fine-scale rel.perm.
			fp.relperm (col_cells.size (col), &wat_sat[0], col_cells[col], &wat_mob[0], 0);
			fp.relperm (col_cells.size (col), &gas_sat[0], col_cells[col], &gas_mob[0], 0);

			// cache the pointers here to avoid indexing in the loop
			double* prm_gas_col = prm_gas[col];
			double* prm_res_col = prm_res[col];
			double* prm_wat_col = prm_wat[col];

			for (int row = 0; row < up.num_rows (col); ++row) {
				// rel.perm. for CO2 when having maximal sat. (only residual brine); this
				// is the rel.perm. for the CO2 that is in the plume
				const double kr_plume = gas_mob[row * num_phases + GAS];

				// rel.perm. of brine, when residual CO2
				const double kr_brine = wat_mob[row + num_phases + WAT];

				// upscaled rel. perm. change for this block; we'll use this to weight
				// the depth fractions when we integrate to get the upscaled rel. perm.
				const double k_factor = lkl[row] / tot_lkl;
				prm_gas_col[row] = k_factor * kr_plume;
				prm_wat_col[row] = k_factor * kr_brine;
				prm_res_col[row] = k_factor * (1 - kr_brine);
			}

			// integrate the derivate to get the upscaled rel. perm.
			up.wgt_dpt (col, prm_gas_col, &prm_gas_int[col][0]);
			up.wgt_dpt (col, prm_wat_col, &prm_wat_int[col][0]);
			up.wgt_dpt (col, prm_res_col, &prm_res_int[col][0]);
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
		// cache this on the outside of the loop
		const int num_phases = numPhases();

		// process each column/cell individually
		for (int i = 0; i < n; ++i) {
			// index (into the upscaled grid) of the column
			const int col = cells[i];

			// get the (upscaled) CO2 saturation
			const double Sg = s[i * num_phases + GAS];

			// get the block number that contains the active interface
			const Elevation intf = intf_elev (col, Sg); // zeta_M

			// rel.perm. for CO2 at this location; simply look up in the
			// table of integrated rel.perm. changes by depth
			const double Krg = up.eval (col, prm_gas_int[col], intf);

			// registered level of maximum CO2 sat. (where there is at least
			// residual CO2
			const Elevation res_lvl = max_gas_elev[col]; // zeta_R

			// rel.perm. for brine at this location; notice that all of
			// our expressions uses the CO2 saturation as parameter
			const double Krw = 1 - (up.eval (col, prm_res_int[col], res_lvl)
			                       +up.eval (col, prm_wat_int[col], intf));

			// assign to output
			kr[i * num_phases + GAS] = Krg;
			kr[i * num_phases + WAT] = Krw;

			// was derivatives requested?
			if (dkrds) {
				// volume available for the mobile liquid/gas: \phi (1-s_{w,r}-s_{g,r})
				const double mob_vol = up.eval (col, mob_mix_vol[col], intf);

				// rel.perm. change for CO2: K^{-1} k_|| k_{r,g}(1-s_{w,r})
				const double prm_chg_gas = up.eval (col, prm_gas[col], intf);

				// possible change in CO2 rel.perm.
				const double dKrg_dSg = upscaled_poro[col] / mob_vol * prm_chg_gas;

				// rel.perm. change for brine: K^{-1} k_|| k_{r,w}(s_{g,r})
				const double prm_chg_wat = up.eval (col, prm_wat[col], intf);

				// possible change in brine rel.perm.
				const double dKrw_dSg = -upscaled_poro[col] / mob_vol * prm_chg_wat;

				// assign to output: since Sw = 1 - Sg, then dkr_g/ds_w = -dkr_g/ds_g
				// viewed as a 2x2 record; the minor index designates the denominator
				// (saturation) and the major index designates the numerator (rel.perm.)
				dkrds[i * (num_phases * num_phases) + num_phases * GAS + GAS] =  dKrg_dSg;
				dkrds[i * (num_phases * num_phases) + num_phases * GAS + WAT] = -dKrg_dSg;
				dkrds[i * (num_phases * num_phases) + num_phases * WAT + GAS] =  dKrw_dSg;
				dkrds[i * (num_phases * num_phases) + num_phases * WAT + WAT] = -dKrw_dSg;
			}
		}
	}

	virtual void capPress (const int n,
	                       const double *s,
	                       const int *cells,
	                       double *pc,
	                       double *dpcds) const {
		// cache this on the outside of the loop; the phase properties
		// are the same in every block
		const int num_phases = numPhases();
		const double dens_gas = density ()[GAS];
		const double dens_wat = density ()[WAT];
		const double dens_diff = dens_gas - dens_wat;

		// wrappers to make sure that we can access this matrix without
		// doing index calculations ourselves
		const rlw_double ts_h (ts.number_of_cells, ts.col_cellpos, ts.h);

		// process each column/cell individually
		for (int i = 0; i < n; ++i) {
			// index (into the upscaled grid) of the column
			const int col = cells[i];

			// get the (upscaled) CO2 saturation
			const double Sg = s[i * num_phases + GAS];

			// get the block number that contains the active interface
			const Elevation intf = intf_elev (col, Sg); // zeta_M

			// heights from top surface to the interface, and to bottom
			const double intf_hgt = up.eval (col, ts_h[col], intf); // \zeta_T - \zeta_M
			const double botm_hgt = ts.h_tot[col]; // \zeta_T - \zeta_B

			// the slopes of the pressure curves are different. the distance
			// between them (at the top for instance) is dependent on where
			// they intersect (i.e. at the interface between the phases). in
			// addition, the brine pressure is measured at the bottom, so we
			// must also add the total height to get down there
			const double hyd_diff = -gravity * (intf_hgt * dens_diff + botm_hgt * dens_wat);

			// total capillary pressure
			const double cap_pres = hyd_diff;

			// assign to output
			pc[i * num_phases + GAS] =  cap_pres;
			pc[i * num_phases + WAT] = -cap_pres;
		}
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
                     const TopSurf& topSurf,
                     const double* grav_vec) {
	// construct real object which contains all the implementation details
	auto_ptr <VertEqProps> props (new VertEqPropsImpl (fineProps,
	                                                   topSurf,
	                                                   grav_vec));

	// client owns pointer to constructed fluid object from this point
	return props.release ();
}
