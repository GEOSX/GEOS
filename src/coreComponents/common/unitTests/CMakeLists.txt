# Specify list of tests
set( gtest_geosx_tests
     testDataTypes.cpp
     testFixedSizeDeque.cpp
     testTypeDispatch.cpp
     testLifoStorage.cpp
     testUnits.cpp )

if ( ENABLE_CALIPER )
  list( APPEND gtest_geosx_tests
        testCaliperSmoke.cpp )
endif()

set( dependencyList ${parallelDeps} common HDF5::HDF5 gtest )

# Add gtest C++ based tests
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )

endforeach()
