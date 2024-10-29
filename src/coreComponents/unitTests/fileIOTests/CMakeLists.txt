# Specify list of tests
set( geosx_fileio_tests testHDFFile.cpp )

set( dependencyList ${parallelDeps} gtest HDF5::HDF5 )

if ( GEOS_BUILD_SHARED_LIBS )
  list( APPEND dependencyList geosx_core )
else()
  list( APPEND dependencyList ${geosx_core_libs} )
endif()

# Add gtest C++ based tests
foreach(test ${geosx_fileio_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )
endforeach()

if ( ENABLE_MPI )

  set(nranks 2)

  set( geosx_fileio_parallel_tests
       testHDFParallelFile.cpp )
  foreach(test ${geosx_fileio_parallel_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name}
                   NUM_MPI_TASKS ${nranks} )
  endforeach()
endif()
