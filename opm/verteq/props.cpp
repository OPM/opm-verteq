#include <opm/verteq/props.hpp>
#include <opm/verteq/topsurf.hpp>
#include <opm/verteq/upscale.hpp>
#include <opm/verteq/utility/exc.hpp>
#include <opm/verteq/utility/runlen.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <algorithm> // fill
#include <cmath> // sqrt
#include <memory> // unique_ptr
#include <vector>
using namespace Opm;
using namespace std;

/**
 * In this module, CO2 is referred to as the "gas" phase even though
 * it is in a supercritical state. This is just to keep an easily
 * identifiable moniker on the variables.
 */
struct VertEqPropsImpl : public VertEqProps {
	/// Get the underlaying fluid information from here
	const IncompPropertiesInterface& fp;

	/// Get the grid information from here
	const TopSurf& ts;

	/// Helper object to do averaging
	const VertEqUpscaler up;

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

	// we assume this ordering of the phases in arrays. these used to be
	// static, but they are now initialized from the phase properties
	const int GAS; // = BlackoilPhases::Liquid;
	const int WAT; // = BlackoilPhases::Aqua;
	static const int NUM_PHASES = 2;
	static const int NUM_PHASES_SQ = NUM_PHASES * NUM_PHASES;

	/// Since Pc_ow = -Pc,wo and dS_o = -dS_w, we can query the capillary
	/// pressure for only the first phase, and then adjust the sign to get
	/// the rest.
	const double phase_sign;

	/// Upscaled porosity; this is \Phi in the papers
	vector <double> upscaled_poro;   // 1/H * int_{\Zeta_B}^{\Zeta_T} \phi dz

	/// Upscaled permeability; this is K in the papers
	vector <double> upscaled_absperm;

	/// Volume fractions of gas phase, used in averaging
	RunLenData <double> res_gas_vol; // \phi S_{n,r}
	RunLenData <double> mob_mix_vol; // \phi (1 - S_{w,r} - S_{n,r})
	RunLenData <double> res_wat_vol; // \phi (1 - S_{w,r})

	/// Volume-of-gas-phase-fraction-weighted depths-fractions
	RunLenData <double> res_gas_dpt; // 1/H * int_{h}^{\zeta_T} \phi S_{n,r} dz
	RunLenData <double> mob_mix_dpt; // 1/H * int_{h}^{\zeta_T} \phi (1 - S_{w,r} - S_{n,r} dz
	RunLenData <double> res_wat_dpt; // 1/H * int_{h}^{\zeta_T} \phi (1 - S_{w,r}) dz

	// we need to keep track of where the plume has been and deposited
	// residual CO2. however, finding the interface is non-trivial and
	// should only be done if we actually see a new maximum of the
	// saturation. this array contains the trigger point for recalc.
	vector <double> max_gas_sat;      // S_{g,max}
	vector <Elevation> max_gas_elev;  // \zeta_R

	virtual void upd_res_sat (const double* snap) {
		// update saturation for each column
		for (int col = 0; col < ts.number_of_cells; ++col) {
			// current CO2 saturation
			const double cur_sat = snap[col * NUM_PHASES + GAS];

			// has it increased? is there more of the plume in this column?
			if (cur_sat > max_gas_sat[col]) {
				// recalculate discretized elevation
				max_gas_elev[col] = res_elev (col, cur_sat);

				// update stored saturation so we test correctly next time
				max_gas_sat[col] = cur_sat;
			}
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
	Elevation res_elev(const int col, const double cur_sat) const {
		// take the highest of the current saturation and the historical
		// highest seen. note that this does NOT update the maximum, allowing
		// this routine to be used for hypothetical saturations
		const double max_sat = std::max (cur_sat, max_gas_sat[col]);

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
		// get the residual interface either by historic values, or if the
		// rel.perm. function is called with a hypothetical new saturation
		// which may be greater.
		const Elevation res_lvl =
		    res_elev (col, std::max (gas_sat, max_gas_sat[col])); // Zeta_R

		// the first term is \Phi * S_g representing the volume of CO2, the
		// second is the integral int_{\zeta_R}^{\zeta_T} \phi s_{g,r} dz,
		// representing the volume of residual CO2; the remainder becomes
		// the mobile CO2 volume. the entire equation is multiplied by 1/H.
		const double res_vol = up.eval (col, res_gas_dpt, res_lvl);
		const double gas_vol = upscaled_poro[col] * gas_sat; // \Phi * S_g
		const double mob_vol = gas_vol - res_vol;

		// lookup to find the height that gives this mobile volume
		const Elevation zeta_M = up.find (col, mob_mix_dpt[col], mob_vol);
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

		// assign which phase is which (e.g. CO2 is first, brine is second)
		// a basic assumption of the vertical equilibrium is that the CO2 is
		// the lightest phase and thus rise to the top of the reservoir
		, GAS (fp.density()[0] < fp.density()[1] ? 0 : 1)
		, WAT (1 - GAS)
		, phase_sign (GAS < WAT ? +1. : -1.)

		// allocate memory for intermediate integrals
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

		// check that we only have two phases
		if (fp.numPhases () != NUM_PHASES) {
			throw OPM_EXC ("Expected %d phases, but got %d", NUM_PHASES, fp.numPhases ());
		}

		// allocate memory to store results for faster lookup later
		upscaled_poro.resize (ts.number_of_cells);
		upscaled_absperm.resize (ts.number_of_cells * PERM_MATRIX_2D);

		// buffers that holds intermediate values for each column;
		// pre-allocate to avoid doing that inside the loop
		vector <double> poro (ts.max_vert_res, 0.); // porosity
		vector <double> kxx (ts.max_vert_res, 0.);  // abs.perm.
		vector <double> kxy (ts.max_vert_res, 0.);
		vector <double> kyy (ts.max_vert_res, 0.);
		vector <double> sgr   (ts.max_vert_res * NUM_PHASES, 0.); // residual CO2
		vector <double> l_swr (ts.max_vert_res * NUM_PHASES, 0.); // 1 - residual brine
		vector <double> lkl (ts.max_vert_res); // magnitude of abs.perm.; k_||

		// saturations and rel.perms. of each phase, assuming maximum filling of...
		vector <double> wat_sat (ts.max_vert_res * NUM_PHASES, 0.); // brine; res. CO2
		vector <double> gas_sat (ts.max_vert_res * NUM_PHASES, 0.); // CO2; res. brine
		vector <double> wat_mob (ts.max_vert_res * NUM_PHASES, 0.); // k_r(S_c=S_{c,r})
		vector <double> gas_mob (ts.max_vert_res * NUM_PHASES, 0.); // k_r(S_c=1-S_{b,r})

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
				const double sgr_ = sgr[row * NUM_PHASES + GAS];
				const double l_swr_ = l_swr[row * NUM_PHASES + GAS];

				// portions of the block that are filled with: residual co2,
				// mobile fluid and residual brine, respectively
				res_gas_col[row] = poro[row] * sgr_;            // \phi*S_{n,r}
				mob_mix_col[row] = poro[row] * (l_swr_ - sgr_); // \phi*(1-S_{w,r}-S_{n_r})
				res_wat_col[row] = poro[row] * l_swr_;          // \phi*(1-S_{w,r}
			}

			// weight the relative depth factor (how close are we towards a
			// completely filled column) with the volume portions. this call
			// to up.wgt_dpt is the same as 1/H int_{h}^{\Zeta_T} ... dz
			up.wgt_dpt (col, &res_gas_col[0], res_gas_dpt);
			up.wgt_dpt (col, &mob_mix_col[0], mob_mix_dpt);
			up.wgt_dpt (col, &res_wat_col[0], res_wat_dpt);

			// now, when we queried the saturation ranges, we got back the min.
			// and max. sat., and when there is min. of one, then there should
			// be max. of the other; however, these data are in different arrays!
			// cross-pick such that we get (min CO2, max brine), (max CO2, min brine)
			// instead of (min CO2, min brine), (max CO2, max brine). this code
			// has no other effect than to satisfy the ordering of items required
			// for the relperm() call
			for (int row = 0; row < col_cells.size (col); ++row) {
				wat_sat[row * NUM_PHASES + GAS] = sgr[row * NUM_PHASES + GAS];
				wat_sat[row * NUM_PHASES + WAT] = l_swr[row * NUM_PHASES + WAT];
				gas_sat[row * NUM_PHASES + GAS] = l_swr[row * NUM_PHASES + GAS];
				gas_sat[row * NUM_PHASES + WAT] = sgr[row * NUM_PHASES + WAT];
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
				const double kr_plume = gas_mob[row * NUM_PHASES + GAS];

				// rel.perm. of brine, when residual CO2
				const double kr_brine = wat_mob[row + NUM_PHASES + WAT];

				// upscaled rel. perm. change for this block; we'll use this to weight
				// the depth fractions when we integrate to get the upscaled rel. perm.
				const double k_factor = lkl[row] / tot_lkl;
				prm_gas_col[row] = k_factor * kr_plume;
				prm_wat_col[row] = k_factor * kr_brine;
				prm_res_col[row] = k_factor * (1 - kr_brine);
			}

			// integrate the derivate to get the upscaled rel. perm.
			up.wgt_dpt (col, prm_gas_col, prm_gas_int);
			up.wgt_dpt (col, prm_wat_col, prm_wat_int);
			up.wgt_dpt (col, prm_res_col, prm_res_int);
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
		return NUM_PHASES;
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
		// process each column/cell individually
		for (int i = 0; i < n; ++i) {
			// index (into the upscaled grid) of the column
			const int col = cells[i];

			// get the (upscaled) CO2 saturation
			const double Sg = s[i * NUM_PHASES + GAS];

			// get the block number that contains the active interface
			const Elevation intf = intf_elev (col, Sg); // zeta_M

			// rel.perm. for CO2 at this location; simply look up in the
			// table of integrated rel.perm. changes by depth
			const double Krg = up.eval (col, prm_gas_int, intf);

			// registered level of maximum CO2 sat. (where there is at least
			// residual CO2
			const Elevation res_lvl =
			    res_elev (col, std::max (max_gas_sat[col], Sg)); // zeta_R

			// rel.perm. for brine at this location; notice that all of
			// our expressions uses the CO2 saturation as parameter
			const double Krw = 1 - (up.eval (col, prm_res_int, res_lvl)
			                       +up.eval (col, prm_wat_int, intf));

			// assign to output
			kr[i * NUM_PHASES + GAS] = Krg;
			kr[i * NUM_PHASES + WAT] = Krw;

			// was derivatives requested?
			if (dkrds) {
				// volume available for the mobile liquid/gas: \phi (1-s_{w,r}-s_{g,r})
				const double mob_vol = up.eval (col, mob_mix_vol, intf);

				// rel.perm. change for CO2: K^{-1} k_|| k_{r,g}(1-s_{w,r})
				const double prm_chg_gas = up.eval (col, prm_gas, intf);

				// possible change in CO2 rel.perm.
				const double dKrg_dSg = upscaled_poro[col] / mob_vol * prm_chg_gas;

				// rel.perm. change for brine: K^{-1} k_|| k_{r,w}(s_{g,r})
				const double prm_chg_wat = up.eval (col, prm_wat, intf);

				// possible change in brine rel.perm.
				const double dKrw_dSg = -upscaled_poro[col] / mob_vol * prm_chg_wat;

				// assign to output: since Sw = 1 - Sg, then dkr_g/ds_w = -dkr_g/ds_g
				// viewed as a 2x2 record; the minor index designates the denominator
				// (saturation) and the major index designates the numerator (rel.perm.)
				dkrds[i * NUM_PHASES_SQ + NUM_PHASES * GAS + GAS] =  dKrg_dSg;
				dkrds[i * NUM_PHASES_SQ + NUM_PHASES * GAS + WAT] = -dKrg_dSg;
				dkrds[i * NUM_PHASES_SQ + NUM_PHASES * WAT + GAS] =  dKrw_dSg;
				dkrds[i * NUM_PHASES_SQ + NUM_PHASES * WAT + WAT] = -dKrw_dSg;
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
		const double dens_gas = density ()[GAS];
		const double dens_wat = density ()[WAT];
		const double dens_diff = dens_gas - dens_wat;

		// wrappers to make sure that we can access this matrix without
		// doing index calculations ourselves
		const rlw_double ts_h (ts.number_of_cells, ts.col_cellpos, ts.h);
		const rlw_double ts_dz (ts.number_of_cells, ts.col_cellpos, ts.dz);
		const rlw_int col_cells (ts.number_of_cells, ts.col_cellpos, ts.col_cells);

		// process each column/cell individually
		for (int i = 0; i < n; ++i) {
			// index (into the upscaled grid) of the column
			const int col = cells[i];

			// get the (upscaled) CO2 saturation
			const double Sg = s[i * NUM_PHASES + GAS];

			// get the block number that contains the active interface
			const Elevation intf = intf_elev (col, Sg); // zeta_M

			// heights from top surface to the interface, and to bottom
			const double intf_hgt = up.eval (col, ts_h, intf); // \zeta_T - \zeta_M

			// the slopes of the pressure curves are different. the distance
			// between them (at the top for instance) is dependent on where
			// they intersect (i.e. at the interface between the phases). if
			// the coordinate system is tilted, we assume that the 'gravity'
			// scalar here is the inner product between the vertical axis and
			// the real gravity vector.
			const double hyd_diff = -gravity * (intf_hgt * dens_diff);

			// find the fine-scale element that holds the interface; we already
			// know the relative index in the column; ask the top surface for
			// global identity
			const int glob_id = col_cells[col][intf.block()];

			// find the entry pressure in this block. this code could
			// be optimized so it only called the capillary pressure
			// function for the fine-scale properties once instead of
			// inside the loop, but that would require us to allocate
			// arrays to hold all input and output, instead of just using
			// local variables. BTW; why the number of outputs?
			double fine_sat[NUM_PHASES];
			double fine_pc[NUM_PHASES];               // entry pressures
			double fine_dpc[NUM_PHASES_SQ];           // derivatives
			fine_sat[GAS] = intf.fraction ();
			fine_sat[WAT] = 1 - fine_sat[GAS];
			fp.capPress (1, fine_sat, &glob_id, fine_pc, fine_dpc);

			// total capillary pressure. the fine scale entry pressure is
			// a wedge between the slopes of the hydrostatic pressures.
			const double fine_pc_GAS = phase_sign * fine_pc[0];
			const double cap_pres = fine_pc_GAS + hyd_diff;

			// assign to output; only the first phase is set, the other should
			// be set to zero (?), see method SimpleFluid2pWrappingProps::pc in
			// opm/core/transport/implicit/SimpleFluid2pWrappingProps_impl.hpp
			pc[i * NUM_PHASES + 0] = phase_sign * cap_pres;
			pc[i * NUM_PHASES + 1] = 0.;

			// interested in the derivatives of the capillary pressure as well?
			if (dpcds) {
				// volume available for the mobile liquid/gas: \phi (1-s_{w,r}-s_{g,r})
				const double mob_vol = up.eval (col, mob_mix_vol, intf);

				// change of interface height per of upscaled saturation; d\zeta_M/dS
				const double dh_dSg = -(ts.h_tot[col] * upscaled_poro[col]) / mob_vol;

				// change of hydrostatic pressure diff per change in interface height
				const double hyd_dPc_dh = -gravity * dens_diff; // dPc/d\zeta_M

				// change in entry pressure per *fine* saturation; notice that only one
				// of the derivatives is set; see the code below for dpcds for the sign
				const double dpe_dsg = GAS < WAT ?
				      +fine_dpc[NUM_PHASES * GAS + GAS] :
				      -fine_dpc[NUM_PHASES * WAT + WAT] ;

				// change in fine saturation per interface height (in this block)
				const double dsg_dh = 1 / ts_dz[col][intf.block()];

				// derivative with respect to upscaled saturation
				const double dPc_dSg = (dpe_dsg * dsg_dh + hyd_dPc_dh) * dh_dSg;

				// assign to output: since Sw = 1 - Sg, then dpc_g/ds_w = -dkr_g/ds_g
				// viewed as a 2x2 record; the minor index designates the denominator
				// (saturation) and the major index designates the numerator (rel.perm.)
				// here too (like for pc) only the first phase is set, the others should
				// have the magic value zero hard-coded (?)
				dpcds[i * NUM_PHASES_SQ + NUM_PHASES * 0 + 0] = phase_sign * dPc_dSg;
				dpcds[i * NUM_PHASES_SQ + NUM_PHASES * 0 + 1] = 0.;
				dpcds[i * NUM_PHASES_SQ + NUM_PHASES * 1 + 0] = 0.;
				dpcds[i * NUM_PHASES_SQ + NUM_PHASES * 1 + 1] = 0.;
			}
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

		// this is just to avoid warnings about unused variable
		static_cast <void> (cells);
	}

	virtual void upscale_pressure (const double* coarseSaturation,
	                               const double* finePressure,
	                               double* coarsePressure) {
		// pressure locations we'll have to relate to
		static const int    FIRST_BLOCK = 0;    // relative index in the column
		static const double HALFWAY     = 0.5;  // center of the block

		// incompressible means that the density is the same everywhere
		// we can thus cache the phase properties outside of the loop
		const double gas_dens = density ()[GAS];
		const double wat_dens = density ()[WAT];

		// helper object to get the index (into the pressure array) and
		// the height of elements in a column
		const rlw_double ts_dz (ts.number_of_cells, ts.col_cellpos, ts.dz);
		const rlw_int col_cells (ts.number_of_cells, ts.col_cellpos, ts.col_cells);

		// upscale each column separately. assume that something like the
		// EQUIL keyword has been used in the Eclipse file and that the
		// pressures are already in equilibrium. thus, we only need to
		// extract the pressure at the reference point (top surface)
		for (int col = 0; col < col_cells.cols (); ++col) {
			// location of the brine-co2 phase contact
			const double gas_sat = coarseSaturation[col * NUM_PHASES + GAS];
			const Elevation& intf_lvl = intf_elev (col, gas_sat);

			// what fraction of the first block from the pressure point (halfway)
			// up to the top surface is of each of the phases? if the interface
			// is below the first block, or if it is further down than halfway,
			// then everything, otherwise the fraction (less than 0.5)
			double gas_frac = (intf_lvl.block () == FIRST_BLOCK) ?
				min (HALFWAY, intf_lvl.fraction ()) : HALFWAY;

			// id of the upper-most block of this column. if there is no
			// blocks, then the TopSurf object wouldn't generate a column.
			const int block_id = col_cells[col][FIRST_BLOCK];

			// height of the uppermost block (twice the distance from the top
			// to the center
			const double hgt = ts_dz[col][FIRST_BLOCK];

			// get the pressure in the middle of this block
			const double mid_pres = finePressure[block_id];

			// pressure at the reference point; adjust hydrostatically for
			// those phases that are on the way up from the center of the first
			// block in the column.
			const double ref_pres = mid_pres - gravity * hgt *
			    (gas_frac * gas_dens + (HALFWAY - gas_frac) * wat_dens);

			// Eclipse uses non-aquous pressure (see Variable Sets in Formulation
			// of the Equations in the Technical Description) as the main unknown
			// in the pressure equation; there is assumed continuity at the
			// contact, so the pressure at the top should always be a CO2 pressure
			coarsePressure[col] = ref_pres;
		}
	}

	virtual void upscale_saturation (const double* fineSaturation,
	                                 double* coarseSaturation) {
		// pointer to all porosities in the fine grid
		const double* fine_poro = fp.porosity ();

		// allocate memory outside of the loop
		vector <double> phi (ts.max_vert_res, 0.); // fine porosity
		vector <double> sg  (ts.max_vert_res, 0.); // fine saturation
		vector <double> pvg (ts.max_vert_res, 0.); // fine pore volume

		// use this object to find the actual number of columns
		rlw_int colcellpos (ts.number_of_cells, ts.col_cellpos, ts.col_cells);

		// upscale column by column
		for (int col = 0; col < ts.number_of_cells; ++col) {
			// porosities for the column. this is the same code as
			// early in the constructor, but we don't save all this
			// data because we don't need it very often.
			up.gather (col, &phi[0], fine_poro, 1, 0);

			// pick out the CO2 saturation from the initial values
			up.gather (col, &sg[0], fineSaturation, NUM_PHASES, GAS);

			// pore-volume occupied by that one phase is the product of
			// volume, porosity and saturation; volume/height is part of
			// the upscaling operator - we assume that every part of the
			// column has the same area (violation of this assumption may
			// cause slight mass balance problems)
			for (int i = 0; i < colcellpos.size (col); ++i) {
				pvg[i] = phi[i] * sg[i];
			}

			// get the sum and update output. notice that we only need
			// the total amount of CO2 in the column
			const double col_porevol = up.dpt_avg (col, &pvg[0]);
			const double upscaled_Sg = col_porevol / upscaled_poro[col];
			coarseSaturation[col * NUM_PHASES + GAS] = upscaled_Sg;
			coarseSaturation[col * NUM_PHASES + WAT] = 1 - upscaled_Sg;
		}
	}

	virtual void downscale_saturation (const double* coarseSaturation,
	                                   double* fineSaturation) {
		// scratch vectors that will hold the minimum and maximum, resp.
		// CO2 saturation. we could get these from res_xxx_vol, but then
		// we would have to dig up the porosity for each cell and divide
		// which is not necessarily faster. if we save this data in the
		// object itself, it may add to memory pressure; I assume instead
		// that it is not expensive for the underlaying properties object
		// to deliver these values on-demand.
		vector <double> sgr   (ts.max_vert_res * NUM_PHASES, 0.); // residual CO2
		vector <double> l_swr (ts.max_vert_res * NUM_PHASES, 0.); // 1 - residual brine

		// indexing object that helps us find the cell in a particular column
		const rlw_int col_cells (ts.number_of_cells, ts.col_cellpos, ts.col_cells);

		// downscale each column individually
		for (int col = 0; col < ts.number_of_cells; ++col) {
			// current height of mobile CO2
			const double gas_hgt = coarseSaturation[col * NUM_PHASES + GAS];

			// make sure residual area of CO2 is up to speed; this ensures
			// that we're looking at current data in the max_gas_elev member
			if (gas_hgt > max_gas_sat[col]) {
				throw OPM_EXC ("Call upd_res_sat before downscale_saturation");
			}

			// height of the interface of residual and mobile CO2, resp.
			const Elevation res_gas = max_gas_elev[col];         // zeta_R
			const Elevation mob_gas = intf_elev (col, gas_hgt);  // zeta_M

			// query the fine properties for the residual saturations; notice
			// that only every other item holds the value for CO2
			const int* ids = col_cells[col];
			fp.satRange (col_cells.size (col), ids, &sgr[0], &l_swr[0]);

			// fill the number of whole blocks which contain mobile CO2 and
			// only residual water (maximum CO2)
			for (int row = 0; row < mob_gas.block (); ++row) {
				const double gas_sat = l_swr[row * NUM_PHASES + GAS];
				const int block = ids[row];
				fineSaturation[block * NUM_PHASES + GAS] = gas_sat;
				fineSaturation[block * NUM_PHASES + WAT] = 1 - gas_sat;
			}

			// then fill the number of *whole* blocks which contain only
			// residual CO2. we start out in the block that was not filled
			// with mobile CO2, i.e. these only fill the *extra* blocks
			// where the plume once was but is not anymore
			for (int row = mob_gas.block(); row < res_gas.block(); ++row) {
				const double gas_sat = sgr[row * NUM_PHASES + GAS];
				const int block = ids[row];
				fineSaturation[block * NUM_PHASES + GAS] = gas_sat;
				fineSaturation[block * NUM_PHASES + WAT] = 1 - gas_sat;
			}

			// fill the remaining of the blocks in the column with pure brine
			for (int row = res_gas.block(); row < col_cells.size (col); ++row) {
				const int block = ids[row];
				fineSaturation[block * NUM_PHASES + GAS] = 0.;
				fineSaturation[block * NUM_PHASES + WAT] = 1.;
			}

			// adjust the block with the mobile/residual interface with its
			// fraction of mobile CO2. since we only have a resolution of one
			// block this sharp interface will only be seen on the visualization
			// as a slightly differently colored block. only do this if there
			// actually is a partially filled block.
			const int intf_block = ids[mob_gas.block ()];
			if (intf_block != col_cells.size(col)) {
				// there will already be residual gas in this block thanks to the
				// loop above; we must only fill a fraction of it with mobile gas,
				// which is the difference between the maximum and minimum filling
				const double intf_gas_sat_incr = mob_gas.fraction () *
				    (l_swr[intf_block * NUM_PHASES + GAS]
				    - sgr[intf_block * NUM_PHASES + GAS]);
				fineSaturation[intf_block * NUM_PHASES + GAS] += intf_gas_sat_incr;
				// we could have written at the brine saturations afterwards to
				// avoid this extra adjustment, but the data locality will be bad
				fineSaturation[intf_block * NUM_PHASES + WAT] -= intf_gas_sat_incr;
			}

			// do the same drill, but with the fraction of where the residual
			// zone ends (the outermost historical edge of the plume)
			const int res_block = ids[res_gas.block()];
			if (res_block != col_cells.size(col)) {
				const double res_gas_sat_incr = res_gas.fraction() *
					  sgr[res_block * NUM_PHASES + GAS];
				fineSaturation[res_block * NUM_PHASES + GAS] += res_gas_sat_incr;
				fineSaturation[res_block * NUM_PHASES + WAT] -= res_gas_sat_incr;
			}
		}
	}

	virtual void downscale_pressure (const double* coarseSaturation,
	                                 const double* coarsePressure,
	                                 double* finePressure) {
		// pressure locations we'll have to relate to
		static const double HALFWAY     = 0.5;  // center of the block

		// helper object to get the index (into the pressure array) and
		// the height of elements in a column
		const rlw_double ts_h (ts.number_of_cells, ts.col_cellpos, ts.h);
		const rlw_double ts_dz (ts.number_of_cells, ts.col_cellpos, ts.dz);
		const rlw_int col_cells (ts.number_of_cells, ts.col_cellpos, ts.col_cells);

		// incompressible means that the density is the same everywhere
		// we can thus cache the phase properties outside of the loop
		const double gas_dens = density ()[GAS];
		const double wat_dens = density ()[WAT];

		for (int col = 0; col < col_cells.cols (); ++col) {
			// location of the brine-co2 phase contact
			const double gas_sat = coarseSaturation[col * NUM_PHASES + GAS];
			const Elevation& intf_lvl = intf_elev (col, gas_sat);

			// get the pressure difference between the phases at top of this column
			const double sat[NUM_PHASES] = { gas_sat, 1-gas_sat };
			double pres_diff;
			capPress (1, sat, &col, &pres_diff, 0);

			// get the reference phase pressure at the top; notice that the CO2
			// pressure is the largest so we subtract the difference
			const double gas_ref = coarsePressure[col];
			const double wat_ref = gas_ref - pres_diff;

			// are we going to include the block with the interface
			const int incl_intf = intf_lvl.fraction () >= HALFWAY ? 1 : 0;
			const int num_gas_rows = intf_lvl.block () + incl_intf;

			// write all CO2 pressure blocks
			for (int row = 0; row < num_gas_rows; ++row) {
				// height of block center
				const double hgt = ts_h[col][row] + HALFWAY * ts_dz[col][row];

				// hydrostatically get the pressure for this block
				const double gas_pres = gas_ref + gravity * hgt * gas_dens;
				const int block = col_cells[col][row];

				// (scatter) write to output array
				finePressure[block] = gas_pres;
			}

			// then write the brine blocks, starting where we left off
			for (int row = num_gas_rows; row < col_cells.size (col); ++row) {
				// height of block center
				const double hgt = ts_h[col][row] + HALFWAY * ts_dz[col][row];

				// hydrostatically get the pressure for this block
				const double wat_pres = wat_ref + gravity * hgt * wat_dens;
				const int block = col_cells[col][row];

				// (scatter) write to output array
				finePressure[block] = wat_pres;
			}
		}
	}
};

VertEqProps*
VertEqProps::create (const IncompPropertiesInterface& fineProps,
                     const TopSurf& topSurf,
                     const double* grav_vec) {
	// construct real object which contains all the implementation details
	unique_ptr <VertEqProps> props (new VertEqPropsImpl (fineProps,
	                                                     topSurf,
	                                                     grav_vec));

	// client owns pointer to constructed fluid object from this point
	return props.release ();
}
