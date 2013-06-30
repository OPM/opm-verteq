// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0
#include <opm/verteq/nav.hpp>
#include <opm/verteq/props.hpp>
#include <opm/verteq/topsurf.hpp>
#include <opm/verteq/verteq.hpp>
#include <opm/verteq/utility/exc.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/wells.h>
#include <memory>           // auto_ptr

using namespace Opm;
using namespace Opm::parameter;
using namespace std;

// Actual implementation of the upscaling
struct VertEqImpl : public VertEq {
	// this pointer needs special handling to dispose; use a zero pointer
	// to signal that it has not been initialized properly (probably some
	// other component which threw an exception)
	Wells* w;
	VertEqImpl () : w (0) {}
	virtual ~VertEqImpl () {
		if (w) {
			destroy_wells (w);
		}
	}
	void init (const UnstructuredGrid& fullGrid,
	           const IncompPropertiesInterface& fullProps,
	           const Wells* wells,
	           const double* gravity);
	// public methods defined in the interface
	virtual const UnstructuredGrid& grid();
	virtual const Wells* wells();
	virtual const IncompPropertiesInterface& props();
	virtual void upscale (const TwophaseState& fineScale,
	                      TwophaseState& coarseScale);
	virtual void notify (const TwophaseState& coarseScale);

	auto_ptr <TopSurf> ts;
	auto_ptr <VertEqProps> pr;
	/**
	 * Translate all the indices in the well list from a full, three-
	 * dimensional grid into the upscaled top surface.
	 */
	void translate_wells ();
};

VertEq*
VertEq::create (const string& title,
                const ParameterGroup& args,
                const UnstructuredGrid& fullGrid,
                const IncompPropertiesInterface& fullProps,
                const Wells* wells,
                const double* gravity) {
	// we don't provide any parameters to do tuning yet
	auto_ptr <VertEqImpl> impl (new VertEqImpl ());
	impl->init (fullGrid, fullProps, wells, gravity);
	return impl.release();
}

void
VertEqImpl::init(const UnstructuredGrid& fullGrid,
                 const IncompPropertiesInterface& fullProps,
                 const Wells* wells,
                 const double* gravity) {
	// generate a two-dimensional upscaling as soon as we get the grid
	ts = auto_ptr <TopSurf> (TopSurf::create (fullGrid));
	pr = auto_ptr <VertEqProps> (VertEqProps::create (fullProps, *ts, gravity));
	// create a separate, but identical, list of wells we can work on
	w = clone_wells(wells);
	translate_wells ();
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

	// TODO: set the initial state from the fine-scale state

	// update the properties from the initial state (the
	// simulation object won't call this method before the
	// first timestep; it assumes that the state is initialized
	// accordingly (which is what we do here now)
	notify (coarseScale);
}

void
VertEqImpl::notify (const TwophaseState& coarseScale) {
	// forward this request to the properties we have stored
	pr->upd_res_sat (&coarseScale.saturation()[0]);
}
