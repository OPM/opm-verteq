OPM Vertical Equilibrium Module
===============================

opm-verteq is an add-on to the simulators in opm-core to do vertical
equilibrium upscaling.

You will need to pass the location of an opm-core build tree or
installation with the `opm-core_ROOT` or `opm_ROOT` parameter to CMake,
or `--with-opm-core` parameter in the Autotools-compatible wrapper.

It is not yet as tested as the other components, and should at this
stage only be used to compare with 3D simulations.

How to wrap on existing simulator
---------------------------------

* You will need to include the header
```c++
#include <opm/verteq/wrapper.hpp>
```

* Wrap the instantiation of the simulator in the wrapper; replace
```c++
Opm::SimulatorIncompTwophase sim (...);
```
with
```c++
Opm::VertEqWrapper <Opm::SimulatorIncompTwophase> sim (...);
```

* The vertical equilibrium wrapper does not write output by itself;
you should register a timestep event handler and write output from it.
The `sync` method must be called in the wrapper to "downscale" the
results to the original grid.

* Pass the header and library directory paths of opm-verteq
to your compiler and linker, respectively.

See the difference between the code in `examples/2d/src/co2_2d.cpp`
and `examples/3d/src/co2_3d.cpp`.
