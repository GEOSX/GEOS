# Specify list of tests
set( LAI_tests
     testDofManager.cpp
     testLAIHelperFunctions.cpp )

set( nranks 2 )

# Add gtest C++ based tests
set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )


foreach(test ${LAI_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${decoratedDependencies} ${tplDependencyList}
                    )

if ( ENABLE_MPI )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} -x ${nranks}
                 NUM_MPI_TASKS ${nranks} )
else()
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endif()
endforeach()
