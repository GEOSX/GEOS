# Specify list of tests
set( testSources
     testGeosxTraits.cpp
     testParsing.cpp
     testUtilities.cpp )

set( dependencyList codingUtilities ${parallelDeps} )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )
# Add gtest C++ based tests
foreach( test ${testSources} )

  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${decoratedDependencies} gtest )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )

endforeach()
