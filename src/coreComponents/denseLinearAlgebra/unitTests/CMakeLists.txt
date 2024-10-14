set( serial_tests
     testBlasLapack.cpp
     testSolveLinearSystem.cpp
     testDenseLASolvers.cpp )

set( dependencyList gtest denseLinearAlgebra )

if ( ENABLE_CUDA )
  list( APPEND dependencyList cuda )
endif()

foreach(test ${serial_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()
