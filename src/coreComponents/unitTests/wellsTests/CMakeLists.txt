# Specify list of tests
set( gtest_geosx_tests
     testReservoirSinglePhaseMSWells.cpp
     testWellEnums.cpp )

set( dependencyList ${parallelDeps} gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  list( APPEND dependencyList geosx_core ${parallelDeps} )
else()
  list( APPEND dependencyList ${geosx_core_libs} ${parallelDeps} )
endif()

if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testReservoirCompositionalMultiphaseMSWells.cpp )

    list( APPEND dependencyList PVTPackage )
endif()

if (TARGET pugixml::pugixml)
  list( APPEND dependencyList pugixml::pugixml )
endif()

if (TARGET pugixml)
  list( APPEND dependencyList pugixml )
endif()

if (TARGET fmt::fmt-header-only)
  list( APPEND dependencyList fmt::fmt-header-only )
endif()

if (TARGET fmt)
  list( APPEND dependencyList fmt )
endif()

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

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
