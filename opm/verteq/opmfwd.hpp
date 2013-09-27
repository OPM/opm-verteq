#ifndef OPM_VERTEQ_OPMFWD_HPP_INCLUDED
#define OPM_VERTEQ_OPMFWD_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

/**
 * Forward declarations necessary to setup a simulation in OPM
 */

struct UnstructuredGrid;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm {

namespace parameter { class ParameterGroup; }
struct Event;
class EventSource;
class IncompPropertiesInterface;
class RockCompressibility;
class WellsManager;
class LinearSolverInterface;
class SimulatorTimer;
class TwophaseState;
class WellsManager;
class WellState;
struct SimulatorReport;

} /* namespace Opm */

#endif /* OPM_VERTEQ_OPMFWD_HPP_INCLUDED */
