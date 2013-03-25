# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap
cmake_minimum_required (VERSION 2.8)
project (co2_2d)
enable_language (CXX)

# build debug by default
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE "Debug")
endif (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# extensions to the build system
# if you use this example to create your own project, you should copy
# the cmake/ folder of an OPM project into the root of your own and then
# remove the double-dots from the path below
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../cmake/Modules")

# always build with debug info
include (UseDebugSymbols)
include (UseWarnings)

# use optimizations if we are doing a release build
include (UseOptimization)

# use tricks to build faster
include (UseFastBuilds)

# this will usually be probed by Findopm-core or specified on the
# command-line. we set it here because this example project is located
# where the heuristics don't find it. if you use this as a template,
# then you should remove this line.
set (opm-core_ROOT "${PROJECT_BINARY_DIR}/../../../opm-core")

include (OpmFind)
find_and_append_package (opm-core REQUIRED)
include_directories (${co2_2d_INCLUDE_DIRS})
link_directories (${co2_2d_LIBRARY_DIRS})
add_definitions (${co2_2d_DEFINITIONS})

# configuration defines necessary for the header files
include (ConfigVars)
configure_vars (
	FILE  CXX  "${PROJECT_BINARY_DIR}/config.h"
	WRITE ${co2_2d_CONFIG_VARS}
	)

# create executables in this directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

# locate source files that makes out program(s)
set (co2_2d_SOURCES "src/co2_2d.cpp")

# setup compilation of this particular executable
add_executable (co2_2d ${co2_2d_SOURCES})
target_link_libraries (co2_2d ${co2_2d_LIBRARIES})

# queue this executable to be stripped
strip_debug_symbols (co2_2d)

# example input file (realize directory names)
configure_file (
  ${PROJECT_SOURCE_DIR}/data/tube.param.in
  ${PROJECT_BINARY_DIR}/bin/tube.param
  )
