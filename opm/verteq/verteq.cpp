// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0
#include <opm/verteq/nav.hpp>
#include <opm/verteq/props.hpp>
#include <opm/verteq/topsurf.hpp>
#include <opm/verteq/upscale.hpp>
#include <opm/verteq/verteq.hpp>
#include <opm/verteq/utility/exc.hpp>
#include <opm/core/pressure/flow_bc.h>
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif /* __clang__ */
#include <opm/core/simulator/initState.hpp>
#ifdef __clang__
#pragma clang diagnostic pop
#endif /* __clang__ */
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/wells.h>
#include <memory>           // unique_ptr

using namespace Opm;
using namespace Opm::parameter;
using namespace std;

// Actual implementation of the upscaling
struct VertEqImpl : public VertEq {
	// this pointer needs special handling to dispose; use a zero pointer
	// to signal that it has not been initialized properly (probably some
	// other component which threw an exception)
	Wells* w;
	FlowBoundaryConditions* bnd_cond;
	VertEqImpl () : w (0), bnd_cond (0) {}
	virtual ~VertEqImpl () {
		if (w) {
			destroy_wells (w);
		}
		if (bnd_cond) {
			flow_conditions_destroy (bnd_cond);
		}
	}
	void init (const UnstructuredGrid& fullGrid,
	           const IncompPropertiesInterface& fullProps,
	           const Wells* wells,
	           const vector<double>& fullSrc,
	           const FlowBoundaryConditions* fullBcs,
	           const double* fullGravity);
	// public methods defined in the interface
	virtual const UnstructuredGrid& grid();
	virtual const Wells* wells();
	virtual const IncompPropertiesInterface& props();
	virtual void upscale (const TwophaseState& fineScale,
	                      TwophaseState& coarseScale);
	virtual void downscale (const TwophaseState &coarseScale,
	                        TwophaseState &fineScale);
	virtual void notify (const TwophaseState& coarseScale);

	unique_ptr <TopSurf> ts;
	unique_ptr <VertEqProps> pr;
	/**
	 * Translate all the indices in the well list from a full, three-
	 * dimensional grid into the upscaled top surface.
	 */
	void translate_wells ();

	// source terms in the upscaled grid
	vector<double> coarseSrc;
	void sum_sources (const vector<double>& fullSrc);
	virtual const vector<double>& src ();

	// boundary conditions
	void assert_noflow (const FlowBoundaryConditions* bcs);
	virtual const FlowBoundaryConditions* bcs ();

	// gravity
	const double* grav_vec;
	virtual const double* gravity ();
};

VertEq*
VertEq::create (const string& title,
                const ParameterGroup& args,
                const UnstructuredGrid& fullGrid,
                const IncompPropertiesInterface& fullProps,
                const Wells* wells,
                const vector<double>& fullSrc,
                const FlowBoundaryConditions* fullBcs,
                const double* fullGravity) {
	// this is just to avoid warnings about unused variables
	static_cast <void> (title);
	static_cast <void> (args);

	// we don't provide any parameters to do tuning yet
	unique_ptr <VertEqImpl> impl (new VertEqImpl ());
	impl->init (fullGrid, fullProps, wells, fullSrc, fullBcs, fullGravity);
	return impl.release();
}

void
VertEqImpl::init(const UnstructuredGrid& fullGrid,
                 const IncompPropertiesInterface& fullProps,
                 const Wells* wells,
                 const vector<double>& fullSrc,
                 const FlowBoundaryConditions* fullBcs,
                 const double* fullGravity) {
	// store a pointer to the original gravity vector passed to us
	grav_vec = fullGravity;

	// generate a two-dimensional upscaling as soon as we get the grid
	ts = unique_ptr <TopSurf> (TopSurf::create (fullGrid));
	pr = unique_ptr <VertEqProps> (VertEqProps::create (fullProps, *ts, grav_vec));
	// create a separate, but identical, list of wells we can work on
	w = clone_wells(wells);
	translate_wells ();
	// sum the volumetric sources in each column
	sum_sources (fullSrc);
	// verify that we haven't specified anything than no-flow boundary
	// conditions (these are the only we support currently)
	// TODO: This should be replaced with code that reads through the
	//       grid and map to proper 2D boundary conditions
	assert_noflow (fullBcs);
	// rely on the fact that no boundary conditions means no-flow
	bnd_cond = flow_conditions_construct (0);
}

const double*
VertEqImpl::gravity () {
	// simply return the original two first items; the underlaying simulator
	// cannot "see" the last dimension, because it only know of two elements.
	return grav_vec;
}

const FlowBoundaryConditions*
VertEqImpl::bcs () {
	// return a set of predefined no-flow conditions
	return bnd_cond;
}

static const char* bc_names[] = {
  "no-flow",
  "pressure",
  "flux"
};

void
VertEqImpl::assert_noflow (const FlowBoundaryConditions* bcs) {
	// loop through the number of boundary conditions specified,
	// checking each individually
	for (size_t i = 0; i < bcs->nbc; ++i) {
		if (bcs->type[i] != BC_NOFLOW) {
			// there is no (portable) format for size_t
			const unsigned long ndx = static_cast <unsigned long> (i);
			throw OPM_EXC ("Boundary condition %lu is %s, only no-flow supported",
			               ndx, bc_names[bcs->type[i]]);
		}
	}
}

void
VertEqImpl::sum_sources (const vector<double>& fineSrc) {
	// helper object that does most of the heavy lifting
	VertEqUpscaler up (*ts);

	// there should be one source term (which may be zero) for each block
	// in the upscaled grid. this may not already be allocated, so ensure
	// that it has the correct size before we begin.
	coarseSrc.resize (ts->number_of_cells, 0.);

	// upscale the source term in each column. since the source term is
	// a volumetric flux, it is a simple addition to get that for the entire
	// column -- no weighting is needed
	for (int col = 0; col < ts->number_of_cells; ++col) {
		coarseSrc[col] = up.sum (col, &fineSrc[0]);
	}
}

const vector<double>&
VertEqImpl::src () {
	return this->coarseSrc;
}

void
VertEqImpl::translate_wells () {
	// number of perforations; we assume that each well is specified with
	// only one perforation, so that there is no column which will end up
	// with more than one well.
	const int num_perfs = w->well_connpos[w->number_of_wells];

	// these flags are used to see if we already have a well in each column
	// a more advanced implementation could perhaps join wells if appropriate.
	vector <int> perforated (ts->number_of_cells, Cart2D::NO_ELEM);

	// translate the index of each well
	for (int i = 0; i < num_perfs; ++i) {
		// three-dimensional placement of the well
		const int fine_id = w->well_cells[i];

		// corresponding position in the two-dimensional grid
		const int coarse_id = ts->fine_col[fine_id];

		// sanity check: do we already have a well here? otherwise mark the
		// spot as taken.
		if (perforated[coarse_id] != Cart2D::NO_ELEM) {
			throw OPM_EXC ("Error translating well %d; column %d is already "
			               "perforated with well in fine cell with id %d",
			               i, coarse_id, fine_id);
		}
		else {
			perforated[coarse_id] = fine_id;
		}

		// overwrite (!) the cell identifier with the 2D one; the list
		// is gradually turned into an upscaled version
		w->well_cells[i] = coarse_id;

		// TODO: Well productivity index is dependent on the drawdown, and
		// the drawdown is dependent on the surrounding reservoir pressure
		// which will change when we are using the upscaled version (it is
		// now at the bottom and not in the height of the well). should the
		// well productivity index be adjusted?
	}
}

const UnstructuredGrid&
VertEqImpl::grid () {
	// simply return the standard part of the grid
	return *(ts.get ());
}

const Wells*
VertEqImpl::wells () {
	// simply return our own list of wells we have translated
	return w;
}

const IncompPropertiesInterface&
VertEqImpl::props () {
	// simply return the standard part of the grid
	return *(pr.get ());
}

void
VertEqImpl::upscale (const TwophaseState& fineScale,
                     TwophaseState& coarseScale) {
	// dimension state object to the top grid
	coarseScale.init (*ts, pr->numPhases ());

	// upscale pressure and saturation to find the initial state of
	// the two-dimensional domain. we only need to set the pressure
	// and saturation, the flux is an output field. these methods
	// are handled by the props class, since it already has access to
	// the densities and weights.
	pr->upscale_saturation (&fineScale.saturation ()[0],
	                        &coarseScale.saturation ()[0]);
	pr->upd_res_sat (&coarseScale.saturation ()[0]);
	pr->upscale_pressure (&coarseScale.saturation ()[0],
	                      &fineScale.pressure ()[0],
	                      &coarseScale.pressure ()[0]);

	// use the regular helper method to initialize the face pressure
	// since it is implemented in the header, we have access to it
	// even though it is in an anonymous namespace!
	initFacePressure (this->grid(), coarseScale);

	// update the properties from the initial state (the
	// simulation object won't call this method before the
	// first timestep; it assumes that the state is initialized
	// accordingly (which is what we do here now)
	notify (coarseScale);
}

void
VertEqImpl::downscale (const TwophaseState &coarseScale,
                       TwophaseState &fineScale) {
	// assume that the fineScale storage is already initialized
	if (!fineScale.pressure().size() == ts->number_of_cells) {
		throw OPM_EXC ("Fine scale state is not dimensioned correctly");
	}

	// properties object handle the actual downscaling since it
	// already has the information about the interface.
	// update the coarse saturation *before* we downscale to 3D,
	// since we need the residual interface for that.
	pr->upd_res_sat (&coarseScale.saturation ()[0]);
	pr->downscale_saturation (&coarseScale.saturation ()[0],
	                          &fineScale.saturation ()[0]);
	pr->downscale_pressure (&coarseScale.saturation ()[0],
	                        &coarseScale.pressure ()[0],
	                        &fineScale.pressure ()[0]);
}

void
VertEqImpl::notify (const TwophaseState& coarseScale) {
	// forward this request to the properties we have stored
	pr->upd_res_sat (&coarseScale.saturation()[0]);
}
