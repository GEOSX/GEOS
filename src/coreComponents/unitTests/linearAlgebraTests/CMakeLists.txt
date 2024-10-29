# Specify list of tests
set( LAI_tests
     testDofManager.cpp
     testLAIHelperFunctions.cpp )

set( nranks 2 )

# Add gtest C++ based tests
set( dependencyList ${parallelDeps} gtest )

if ( GEOS_BUILD_SHARED_LIBS )
  list( APPEND dependencyList geosx_core HDF5::HDF5 )
else()
  list( APPEND dependencyList ${geosx_core_libs} HDF5::HDF5 )
endif()

if (TARGET pugixml::pugixml)
  set( dependencyList ${dependencyList} pugixml::pugixml )
endif()

if (TARGET pugixml)
  set( dependencyList ${dependencyList} pugixml )
endif()

if (TARGET fmt::fmt-header-only)
  set( dependencyList ${dependencyList} fmt::fmt-header-only )
endif()

if (TARGET fmt)
  set( dependencyList ${dependencyList} fmt )
endif()

foreach(test ${LAI_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList}
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
