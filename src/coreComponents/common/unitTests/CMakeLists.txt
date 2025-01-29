# Specify list of tests
set( gtest_geosx_tests
     testDataTypes.cpp
     testFixedSizeDeque.cpp
     testMpiWrapper.cpp
     testTypeDispatch.cpp
     testLifoStorage.cpp
     testUnits.cpp )

set( gtest_geosx_mpi_tests
     testMpiWrapper.cpp )

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

if( ENABLE_MPI )
  set( nranks 4 )

  foreach( test ${gtest_geosx_mpi_tests} )
    get_filename_component( file_we ${test} NAME_WE )
    set( test_name ${file_we}_mpi )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
  endforeach()
endif()
