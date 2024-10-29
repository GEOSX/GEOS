# Specify list of tests
set( gtest_geosx_tests
     testMimeticInnerProducts.cpp )

set( dependencyList ${parallelDeps} gtest )

if ( GEOS_BUILD_SHARED_LIBS )
  list( APPEND dependencyList geosx_core HDF5::HDF5 )
else()
  list( APPEND dependencyList ${geosx_core_libs} HDF5::HDF5 )
  list( APPEND dependencyList ${geosx_core_libs} HDF5::HDF5 )
endif()

# Add gtest C++ based tests
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    # Guard to prevent GCC (version 8) from giving warnings due
    # to some sort of possible conversion from int to long unsigned.
    # See also discussion here: https://github.com/GEOS-DEV/LvArray/pull/250
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.3.1)
        target_compile_options(${test_name} PRIVATE "-Wno-stringop-overflow")
      endif()
    endif()

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )

endforeach()
