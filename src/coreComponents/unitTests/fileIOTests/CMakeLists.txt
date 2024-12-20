# Specify list of tests
set( geosx_fileio_tests testHDFFile.cpp )

set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

# Add gtest C++ based tests
foreach(test ${geosx_fileio_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

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
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name}
                   NUM_MPI_TASKS ${nranks} )
  endforeach()
endif()
