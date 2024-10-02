# Specify list of tests
set( gtest_geosx_tests
     testReservoirSinglePhaseMSWells.cpp
     testWellEnums.cpp )

set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface )

if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testReservoirCompositionalMultiphaseMSWells.cpp )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )
# Add gtest C++ based tests
foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
