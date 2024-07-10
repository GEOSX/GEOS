# Specify list of tests
set( gtest_geosx_tests
     testWavePropagation.cpp
     testWavePropagationQ2.cpp
     testWavePropagationElasticFirstOrder.cpp
     testWavePropagationDAS.cpp
     testWavePropagationElasticVTI.cpp
     testWavePropagationAttenuation.cpp
     testWavePropagationAcousticFirstOrder.cpp )

set( dependencyList ${parallelDeps} gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  list( APPEND dependencyList geosx_core )
else()
  list( APPEND dependencyList ${geosx_core_libs} )
endif()

if (TARGET pugixml::pugixml)
  list( APPEND dependencyList pugixml::pugixml )
endif()

if (TARGET pugixml)
  list( APPEND dependencyList pugixml )
endif()

if (TARGET fmt::fmt)
  list( APPEND dependencyList fmt::fmt )
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

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
