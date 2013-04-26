// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0
#include <opm/verteq/topsurf.hpp>
#include <opm/verteq/verteq.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <memory>           // auto_ptr

using namespace Opm;
using namespace Opm::parameter;
using namespace std;

// Actual implementation of the upscaling
struct VertEqImpl : public VertEq {
	// public methods defined in the interface
	virtual ~VertEqImpl () {}
	void init (const UnstructuredGrid& fullGrid);
	virtual const UnstructuredGrid& grid();

	auto_ptr <TopSurf> ts;
};

VertEq*
VertEq::create (const string& title,
								const ParameterGroup& args,
								const UnstructuredGrid& fullGrid) {
	// we don't provide any parameters to do tuning yet
	auto_ptr <VertEqImpl> impl (new VertEqImpl ());
	impl->init (fullGrid);
	return impl.release();
}

void
VertEqImpl::init(const UnstructuredGrid& fullGrid) {
	// generate a two-dimensional upscaling as soon as we get the grid
	ts = auto_ptr <TopSurf> (TopSurf::create (fullGrid));
}

const UnstructuredGrid&
VertEqImpl::grid () {
	// simply return the standard part of the grid
	return *(ts.get ());
}
