// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0
#include <opm/verteq/props.hpp>
#include <opm/verteq/topsurf.hpp>
#include <opm/verteq/verteq.hpp>
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
