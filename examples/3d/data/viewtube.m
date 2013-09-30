% inspired by modules/impes/examples/SPE1/deckinputSPE1Example.m
mrstModule add deckformat

% get the data file relative to this script
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'tube.data');

% read grid from the input file
deck  = readEclipseDeck(fn);
deck  = convertDeckUnits(deck);
fluid = initEclipseFluid(deck);
rock  = initEclipseRock(deck);
grid  = initEclipseGrid(deck);
grid  = computeGeometry(grid);

% read wells from the input file
rock  = compressRock(rock, grid.cells.indexMap);
wells = processWells(grid, rock, deck.SCHEDULE.control, 'InnerProduct', 'ip_tpf');

% get initial state from the input file
state = initEclipseState(grid, deck, fluid);
pres  = state.pressure(1);
state.wellSol = initWellSol(wells, pres);

% plot grid and wells
clf;
plotGrid(grid, 'FaceAlpha', .3, 'EdgeAlpha', .1);
plotWell(grid, wells);
view(0, 0);
axis equal tight;
