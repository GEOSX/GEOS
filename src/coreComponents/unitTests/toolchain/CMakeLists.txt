# Specify list of tests
set( gtest_geosx_tests testToolchain.cpp )

# NOTE: we explicitly depend on internal libraries in order to check for the spurious addition of -DNDEBUG flag
set( dependencyList ${parallelDeps}  physicsSolvers )
set( tplDependencyList ${parallelDeps} gtest HDF5::HDF5 optionparser )
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

    # NOTE: we explicitly depend on internal headers in order to check for the spurious addition of -DNDEBUG flag
    target_include_directories( ${test_name} PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
endforeach()
