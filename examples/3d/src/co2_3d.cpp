/* -*- mode: c++; tab-width: 2; indent-tabs-mode: t; truncate-lines: t -*- */
/* vim: set filetype=cpp autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap: */
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/IncompPropertiesFromDeck.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/TwophaseState.hpp>

#include <iostream>

using namespace Opm;
using namespace Opm::parameter;
using namespace std;

int main (int argc, char *argv[]) {
	// read parameters from command-line
	ParameterGroup param (argc, argv, false);

	// parse keywords from the input file specified
	// TODO: requirement that file exists
	const string filename = param.get <string> ("filename");
	cout << "Reading deck: " << filename << endl;
	const EclipseGridParser parser (filename);

	// extract grid from the parse tree
	const GridManager gridMan (parser);
	const UnstructuredGrid& grid = *gridMan.c_grid ();

	// extract fluid, rock and two-phase properties from the parse tree
	IncompPropertiesFromDeck fluid (parser, grid);

	// initial state of the reservoir
	const double gravity [] = { 0., 0., Opm::unit::gravity };
	TwophaseState state;
	initStateFromDeck (grid, fluid, parser, gravity [3], state);

	// if some parameters were unused, it may be that they're spelled wrong
	if (param.anyUnused ()) {
		cerr << "Unused parameters:" << endl;
		param.displayUsage ();
	}

	// done
	return 0;
}
