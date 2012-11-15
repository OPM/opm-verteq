/* -*- mode: c++; tab-width: 2; indent-tabs-mode: t; truncate-lines: t -*- */
/* vim: set filetype=cpp autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap: */
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <iostream>

using namespace Opm::parameter;
using namespace std;

int main (int argc, char *argv[]) {
	// read parameters from command-line
	ParameterGroup param (argc, argv, false);

	// if some parameters were unused, it may be that they're spelled wrong
	if (param.anyUnused ()) {
		cerr << "Unused parameters:" << endl;
		param.displayUsage ();
	}

	// done
	return 0;
}
