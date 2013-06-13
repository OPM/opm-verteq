#include <opm/verteq/props.hpp>
#include <opm/verteq/topsurf.hpp>
#include <opm/verteq/utility/exc.hpp>
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

	VertEqPropsImpl (const IncompPropertiesInterface& fineProps,
	                 const TopSurf& topSurf)
		: fp (fineProps)
		, ts (topSurf) {
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
		throw OPM_EXC ("Not implemented yet");
	}

	virtual const double* permeability () const {
		throw OPM_EXC ("Not implemented yet");
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
		throw OPM_EXC ("Not implemented yet");
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
