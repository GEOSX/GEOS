set( serial_tests
     testLinearSolverParametersEnums.cpp
     testComponentMask.cpp )

set( parallel_tests
     testMatrices.cpp
     testVectors.cpp
     testExternalSolvers.cpp
     testKrylovSolvers.cpp
     testReverseCutHillMcKeeOrdering.cpp )

set( nranks 2 )

set( dependencyList gtest linearAlgebra ${parallelDeps} )

foreach(test ${serial_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()

foreach(test ${parallel_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )
  if ( ENABLE_MPI )
    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
  else()
    geos_add_test( NAME _${test_name}
                   COMMAND ${test_name} )
  endif()
endforeach()
