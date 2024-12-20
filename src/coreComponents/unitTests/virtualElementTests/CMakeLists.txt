# Specify list of tests
set( testSources
     testConformingVirtualElementOrder1.cpp )

set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )


# Add gtest C++ based tests
foreach(test ${testSources})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

    geos_add_test( NAME ${test_name}
                  COMMAND ${test_name} )

    if( ( NOT ${ENABLE_CUDA} ) AND ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" ) )
        set_target_properties(${test_name} PROPERTIES COMPILE_FLAGS "-Wno-everything")
    endif()

endforeach()

if ( ENABLE_CUDA )
#blt_add_target_compile_flags( TO testShapeFunctionKernel
#                              FLAGS "--maxrregcount=64")
endif()
