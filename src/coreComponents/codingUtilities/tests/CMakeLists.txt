# Specify list of tests
set( testSources
     testGeosxTraits.cpp
     testStringUtilities.cpp
     testParsing.cpp
     testUtilities.cpp )

set( dependencyList gtest codingUtilities ${parallelDeps} )

# Add gtest C++ based tests
foreach( test ${testSources} )

  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )

endforeach()
