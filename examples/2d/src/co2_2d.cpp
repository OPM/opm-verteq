/* -*- mode: c++; tab-width: 2; indent-tabs-mode: t; truncate-lines: t -*- */
/* vim: set filetype=cpp autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap: */
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/IncompPropertiesFromDeck.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/core/simulator/SimulatorIncompTwophase.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorOutput.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/verteq/wrapper.hpp>

#include <iostream>
#include <memory>
#include <vector>

using namespace Opm;
using namespace Opm::parameter;
using namespace std;

int main (int argc, char *argv[]) try {
	// read parameters from command-line
	ParameterGroup param (argc, argv, false);

	// parse keywords from the input file specified
	// TODO: requirement that file exists
	const string filename = param.get <string> ("deck_filename");
	cout << "Reading deck: " << filename << endl;
	const Parser deckGenerator;
	auto deck = deckGenerator.parseFile (filename);

	// extract grid from the parse tree
	const GridManager gridMan (deck);
	const UnstructuredGrid& grid = *gridMan.c_grid ();

	// extract fluid, rock and two-phase properties from the parse tree
	IncompPropertiesFromDeck fluid (deck, grid);

	// initial state of the reservoir
	const double gravity [] = { 0., 0., Opm::unit::gravity };
	TwophaseState state;
	initStateFromDeck (grid, fluid, deck, gravity [2], state);

	// setup wells from input, using grid and rock properties read earlier
	auto eclipseState = make_shared <EclipseState> (deck);
	int reportStepIdx = 0;
	WellsManager wells (eclipseState, reportStepIdx, grid, fluid.permeability());
	WellState wellState; wellState.init (wells.c_wells(), state);

	// no sources and no-flow boundary conditions
	vector <double> src (grid.number_of_cells, 0.);
	FlowBCManager bc;

	// run schedule
	SimulatorTimer stepping;
	auto timeMap = make_shared <TimeMap> (deck);
	stepping.init (timeMap);

	// pressure and transport solvers
	LinearSolverFactory linsolver (param);
	VertEqWrapper<SimulatorIncompTwophase> sim (
		param, grid, fluid, 0, wells, src, bc.c_bcs(), linsolver, gravity);

	// write the state at all reporting times
	SimulatorOutput <VertEqWrapper <SimulatorIncompTwophase> > outp (
		param, *deck, *timeMap, grid, stepping, state, wellState, sim); (void) outp;

	// if some parameters were unused, it may be that they're spelled wrong
	if (param.anyUnused ()) {
		cerr << "Unused parameters:" << endl;
		param.displayUsage ();
	}

	// loop solvers until final time has arrived
	sim.run (stepping, state, wellState);

	// done
	return 0;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

