# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up five lists:
#	MAIN_SOURCE_FILES     List of compilation units which will be included in
#	                      the library. If it isn't on this list, it won't be
#	                      part of the library. Please try to keep it sorted to
#	                      maintain sanity.
#
#	TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
#	TEST_DATA_FILES       Files from the source three that should be made
#	                      available in the corresponding location in the build
#	                      tree in order to run tests there.
#
#	EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#	                      build, but which is not part of the library nor is
#	                      run as tests.
#
#	PUBLIC_HEADER_FILES   List of public header files that should be
#	                      distributed together with the library. The source
#	                      files can of course include other files than these;
#	                      you should only add to this list if the *user* of
#	                      the library needs it.

# originally generated with the command:
# find opm -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
	opm/core/eclipse/EclipseGridInspector.cpp
	opm/core/eclipse/EclipseGridParser.cpp
	opm/core/fluid/blackoil/BlackoilPvtProperties.cpp
	opm/core/fluid/BlackoilPropertiesBasic.cpp
	opm/core/fluid/BlackoilPropertiesFromDeck.cpp
	opm/core/fluid/blackoil/SinglePvtDead.cpp
	opm/core/fluid/blackoil/SinglePvtDeadSpline.cpp
	opm/core/fluid/blackoil/SinglePvtInterface.cpp
	opm/core/fluid/blackoil/SinglePvtLiveGas.cpp
	opm/core/fluid/blackoil/SinglePvtLiveOil.cpp
	opm/core/fluid/IncompPropertiesBasic.cpp
	opm/core/fluid/IncompPropertiesFromDeck.cpp
	opm/core/fluid/PvtPropertiesBasic.cpp
	opm/core/fluid/PvtPropertiesIncompFromDeck.cpp
	opm/core/fluid/RockBasic.cpp
	opm/core/fluid/RockCompressibility.cpp
	opm/core/fluid/RockFromDeck.cpp
	opm/core/fluid/SatFuncGwseg.cpp
	opm/core/fluid/SatFuncSimple.cpp
	opm/core/fluid/SatFuncStone2.cpp
	opm/core/fluid/SaturationPropsBasic.cpp
	opm/core/fluid/SaturationPropsFromDeck.cpp
	opm/core/grid.c
	opm/core/grid/cart_grid.c
	opm/core/grid/cornerpoint_grid.c
	opm/core/grid/cpgpreprocess/facetopology.c
	opm/core/grid/cpgpreprocess/geometry.c
	opm/core/grid/cpgpreprocess/mxgrdecl.c
	opm/core/grid/cpgpreprocess/preprocess.c
	opm/core/grid/cpgpreprocess/processgrid.c
	opm/core/grid/cpgpreprocess/uniquepoints.c
	opm/core/GridManager.cpp
	opm/core/linalg/call_umfpack.c
	opm/core/linalg/LinearSolverAGMG.cpp
	opm/core/linalg/LinearSolverFactory.cpp
	opm/core/linalg/LinearSolverInterface.cpp
	opm/core/linalg/LinearSolverIstl.cpp
	opm/core/linalg/LinearSolverUmfpack.cpp
	opm/core/linalg/sparse_sys.c
	opm/core/newwells.c
	opm/core/pressure/cfsh.c
	opm/core/pressure/CompressibleTpfa.cpp
	opm/core/pressure/flow_bc.c
	opm/core/pressure/FlowBCManager.cpp
	opm/core/pressure/fsh.c
	opm/core/pressure/fsh_common_impl.c
	opm/core/pressure/ifsh.c
	opm/core/pressure/IncompTpfa.cpp
	opm/core/pressure/mimetic/hybsys.c
	opm/core/pressure/mimetic/hybsys_global.c
	opm/core/pressure/mimetic/mimetic.c
	opm/core/pressure/msmfem/coarse_conn.c
	opm/core/pressure/msmfem/coarse_sys.c
	opm/core/pressure/msmfem/dfs.c
	opm/core/pressure/msmfem/hash_set.c
	opm/core/pressure/msmfem/ifsh_ms.c
	opm/core/pressure/msmfem/partition.c
	opm/core/pressure/tpfa/cfs_tpfa.c
	opm/core/pressure/tpfa/cfs_tpfa_residual.c
	opm/core/pressure/tpfa/compr_bc.c
	opm/core/pressure/tpfa/compr_quant.c
	opm/core/pressure/tpfa/compr_quant_general.c
	opm/core/pressure/tpfa/compr_source.c
	opm/core/pressure/tpfa/ifs_tpfa.c
	opm/core/pressure/tpfa/trans_tpfa.c
	opm/core/pressure/well.c
	opm/core/simulator/SimulatorCompressibleTwophase.cpp
	opm/core/simulator/SimulatorIncompTwophase.cpp
	opm/core/simulator/SimulatorReport.cpp
	opm/core/simulator/SimulatorTimer.cpp
	opm/core/transport/reorder/DGBasis.cpp
	opm/core/transport/reorder/nlsolvers.c
	opm/core/transport/reorder/reordersequence.cpp
	opm/core/transport/reorder/tarjan.c
	opm/core/transport/reorder/TransportModelCompressibleTwophase.cpp
	opm/core/transport/reorder/TransportModelInterface.cpp
	opm/core/transport/reorder/TransportModelTracerTof.cpp
	opm/core/transport/reorder/TransportModelTracerTofDiscGal.cpp
	opm/core/transport/reorder/TransportModelTwophase.cpp
	opm/core/transport/spu_explicit.c
	opm/core/transport/spu_implicit.c
	opm/core/transport/transport_source.c
	opm/core/utility/miscUtilitiesBlackoil.cpp
	opm/core/utility/miscUtilities.cpp
	opm/core/utility/MonotCubicInterpolator.cpp
	opm/core/utility/parameters/Parameter.cpp
	opm/core/utility/parameters/ParameterGroup.cpp
	opm/core/utility/parameters/ParameterTools.cpp
	opm/core/utility/parameters/ParameterXML.cpp
	opm/core/utility/parameters/tinyxml/tinystr.cpp
	opm/core/utility/parameters/tinyxml/tinyxml.cpp
	opm/core/utility/parameters/tinyxml/tinyxmlerror.cpp
	opm/core/utility/parameters/tinyxml/tinyxmlparser.cpp
	opm/core/utility/parameters/tinyxml/xmltest.cpp
	opm/core/utility/StopWatch.cpp
	opm/core/utility/VelocityInterpolation.cpp
	opm/core/utility/WachspressCoord.cpp
	opm/core/utility/writeECLData.cpp
	opm/core/utility/writeVtkData.cpp
	opm/core/vag_format/vag.cpp
	opm/core/wells/InjectionSpecification.cpp
	opm/core/wells/ProductionSpecification.cpp
	opm/core/wells/WellCollection.cpp
	opm/core/wells/WellsGroup.cpp
	opm/core/wells/WellsManager.cpp
	)

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
	tests/test_dgbasis.cpp
	tests/test_sparsevector.cpp
	tests/test_sparsetable.cpp
	tests/test_velocityinterpolation.cpp
	tests/test_quadratures.cpp
	tests/test_wells.cpp
	tests/test_wachspresscoord.cpp
	tests/test_column_extract.cpp
	tests/test_geom2d.cpp
	tests/test_param.cpp
	)

# originally generated with the command:
# find tests -name '*.xml' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_DATA_FILES
	tests/extratestdata.xml
	tests/testdata.xml
	)

# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
	examples/compute_tof.cpp
	examples/compute_tof_from_files.cpp
	examples/import_rewrite.cpp
	examples/refine_wells.cpp
	examples/scaneclipsedeck.c
	examples/sim_2p_comp_reorder.cpp
	examples/sim_2p_incomp_reorder.cpp
	examples/sim_wateroil.cpp
	examples/spu_2p.cpp
	examples/wells_example.cpp
	tutorials/tutorial1.cpp
	tutorials/tutorial2.cpp
	tutorials/tutorial3.cpp
	tutorials/tutorial4.cpp
	)

# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
	opm/core/transport/ImplicitTransport.hpp
	opm/core/transport/SimpleFluid2pWrapper.hpp
	opm/core/transport/CSRMatrixBlockAssembler.hpp
	opm/core/transport/ImplicitAssembly.hpp
	opm/core/transport/CSRMatrixUmfpackSolver.hpp
	opm/core/transport/GravityColumnSolver_impl.hpp
	opm/core/transport/spu_implicit.h
	opm/core/transport/transport_source.h
	opm/core/transport/SinglePointUpwindTwoPhase.hpp
	opm/core/transport/JacobianSystem.hpp
	opm/core/transport/reorder/TransportModelTwophase.hpp
	opm/core/transport/reorder/TransportModelTracerTof.hpp
	opm/core/transport/reorder/TransportModelCompressibleTwophase.hpp
	opm/core/transport/reorder/DGBasis.hpp
	opm/core/transport/reorder/reordersequence.h
	opm/core/transport/reorder/tarjan.h
	opm/core/transport/reorder/TransportModelInterface.hpp
	opm/core/transport/reorder/nlsolvers.h
	opm/core/transport/reorder/TransportModelTracerTofDiscGal.hpp
	opm/core/transport/GravityColumnSolver.hpp
	opm/core/transport/NormSupport.hpp
	opm/core/transport/spu_explicit.h
	opm/core/eclipse/SpecialEclipseFields.hpp
	opm/core/eclipse/EclipseUnits.hpp
	opm/core/eclipse/EclipseGridParserHelpers.hpp
	opm/core/eclipse/CornerpointChopper.hpp
	opm/core/eclipse/EclipseGridParser.hpp
	opm/core/eclipse/EclipseGridInspector.hpp
	opm/core/GridManager.hpp
	opm/core/linalg/blas_lapack.h
	opm/core/linalg/LinearSolverAGMG.hpp
	opm/core/linalg/LinearSolverFactory.hpp
	opm/core/linalg/LinearSolverUmfpack.hpp
	opm/core/linalg/LinearSolverIstl.hpp
	opm/core/linalg/call_umfpack.h
	opm/core/linalg/sparse_sys.h
	opm/core/linalg/LinearSolverInterface.hpp
	opm/core/newwells.h
	opm/core/vag_format/vag.hpp
	opm/core/doxygen_main.hpp
	opm/core/simulator/SimulatorIncompTwophase.hpp
	opm/core/simulator/WellState.hpp
	opm/core/simulator/SimulatorCompressibleTwophase.hpp
	opm/core/simulator/SimulatorReport.hpp
	opm/core/simulator/SimulatorTimer.hpp
	opm/core/simulator/TwophaseState.hpp
	opm/core/simulator/BlackoilState.hpp
	opm/core/utility/linearInterpolation.hpp
	opm/core/utility/WachspressCoord.hpp
	opm/core/utility/initState_impl.hpp
	opm/core/utility/initState.hpp
	opm/core/utility/NonuniformTableLinear.hpp
	opm/core/utility/Units.hpp
	opm/core/utility/DataMap.hpp
	opm/core/utility/linInt.hpp
	opm/core/utility/parameters/ParameterXML.hpp
	opm/core/utility/parameters/ParameterGroup_impl.hpp
	opm/core/utility/parameters/Parameter.hpp
	opm/core/utility/parameters/ParameterTools.hpp
	opm/core/utility/parameters/ParameterRequirement.hpp
	opm/core/utility/parameters/ParameterMapItem.hpp
	opm/core/utility/parameters/ParameterGroup.hpp
	opm/core/utility/parameters/tinyxml/tinystr.h
	opm/core/utility/parameters/tinyxml/tinyxml.h
	opm/core/utility/parameters/ParameterStrings.hpp
	opm/core/utility/miscUtilitiesBlackoil.hpp
	opm/core/utility/MonotCubicInterpolator.hpp
	opm/core/utility/UniformTableLinear.hpp
	opm/core/utility/writeVtkData.hpp
	opm/core/utility/Average.hpp
	opm/core/utility/Factory.hpp
	opm/core/utility/VelocityInterpolation.hpp
	opm/core/utility/SparseVector.hpp
	opm/core/utility/StopWatch.hpp
	opm/core/utility/buildUniformMonotoneTable.hpp
	opm/core/utility/ErrorMacros.hpp
	opm/core/utility/miscUtilities.hpp
	opm/core/utility/ColumnExtract.hpp
	opm/core/utility/RootFinders.hpp
	opm/core/utility/SparseTable.hpp
	opm/core/utility/have_boost_redef.hpp
	opm/core/utility/writeECLData.hpp
	opm/core/grid/CellQuadrature.hpp
	opm/core/grid/cornerpoint_grid.h
	opm/core/grid/cpgpreprocess/geometry.h
	opm/core/grid/cpgpreprocess/uniquepoints.h
	opm/core/grid/cpgpreprocess/mxgrdecl.h
	opm/core/grid/cpgpreprocess/preprocess.h
	opm/core/grid/cpgpreprocess/grdecl.h
	opm/core/grid/cpgpreprocess/facetopology.h
	opm/core/grid/cart_grid.h
	opm/core/grid/FaceQuadrature.hpp
	opm/core/pressure/TPFAPressureSolver.hpp
	opm/core/pressure/fsh_common_impl.h
	opm/core/pressure/msmfem/dfs.h
	opm/core/pressure/msmfem/coarse_conn.h
	opm/core/pressure/msmfem/hash_set.h
	opm/core/pressure/msmfem/partition.h
	opm/core/pressure/msmfem/ifsh_ms.h
	opm/core/pressure/msmfem/coarse_sys.h
	opm/core/pressure/IncompTpfa.hpp
	opm/core/pressure/HybridPressureSolver.hpp
	opm/core/pressure/fsh.h
	opm/core/pressure/flow_bc.h
	opm/core/pressure/tpfa/compr_quant_general.h
	opm/core/pressure/tpfa/compr_bc.h
	opm/core/pressure/tpfa/trans_tpfa.h
	opm/core/pressure/tpfa/compr_quant.h
	opm/core/pressure/tpfa/cfs_tpfa_residual.h
	opm/core/pressure/tpfa/ifs_tpfa.h
	opm/core/pressure/tpfa/compr_source.h
	opm/core/pressure/tpfa/cfs_tpfa.h
	opm/core/pressure/mimetic/hybsys_global.h
	opm/core/pressure/mimetic/hybsys.h
	opm/core/pressure/mimetic/mimetic.h
	opm/core/pressure/FlowBCManager.hpp
	opm/core/pressure/CompressibleTpfa.hpp
	opm/core/pressure/TPFACompressiblePressureSolver.hpp
	opm/core/wells/InjectionSpecification.hpp
	opm/core/wells/WellsGroup.hpp
	opm/core/wells/WellsManager.hpp
	opm/core/wells/ProductionSpecification.hpp
	opm/core/wells/WellCollection.hpp
	opm/core/fluid/SaturationPropsInterface.hpp
	opm/core/fluid/SaturationPropsFromDeck_impl.hpp
	opm/core/fluid/RockFromDeck.hpp
	opm/core/fluid/BlackoilPropertiesBasic.hpp
	opm/core/fluid/RockCompressibility.hpp
	opm/core/fluid/PvtPropertiesIncompFromDeck.hpp
	opm/core/fluid/IncompPropertiesBasic.hpp
	opm/core/fluid/SaturationPropsFromDeck.hpp
	opm/core/fluid/SaturationPropsBasic.hpp
	opm/core/fluid/BlackoilPropertiesFromDeck.hpp
	opm/core/fluid/BlackoilPropertiesInterface.hpp
	opm/core/fluid/blackoil/SinglePvtInterface.hpp
	opm/core/fluid/blackoil/BlackoilPhases.hpp
	opm/core/fluid/blackoil/SinglePvtDeadSpline.hpp
	opm/core/fluid/blackoil/SinglePvtLiveOil.hpp
	opm/core/fluid/blackoil/SinglePvtConstCompr.hpp
	opm/core/fluid/blackoil/SinglePvtDead.hpp
	opm/core/fluid/blackoil/phaseUsageFromDeck.hpp
	opm/core/fluid/blackoil/BlackoilPvtProperties.hpp
	opm/core/fluid/blackoil/SinglePvtLiveGas.hpp
	opm/core/fluid/SimpleFluid2p.hpp
	opm/core/fluid/IncompPropertiesInterface.hpp
	opm/core/fluid/PvtPropertiesBasic.hpp
	opm/core/fluid/SatFuncSimple.hpp
	opm/core/fluid/SatFuncGwseg.hpp
	opm/core/fluid/IncompPropertiesFromDeck.hpp
	opm/core/fluid/RockBasic.hpp
	opm/core/fluid/SatFuncStone2.hpp
	opm/core/grid.h
	opm/core/well.h
	opm/core/GridAdapter.hpp
	)
