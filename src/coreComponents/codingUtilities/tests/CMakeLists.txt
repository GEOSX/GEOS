# Specify list of tests
set( testSources
     testGeosxTraits.cpp
     testParsing.cpp
     testRTTypes.cpp
     testUtilities.cpp )

set(mock_object_headers
    BaseClass.hpp
    DerivedFinalClass.hpp
)

set(mock_object_impl
    BaseClass.cpp
    DerivedFinalClass.cpp
)

set( dependencyList codingUtilities ${parallelDeps})

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )
# Add gtest C++ based tests
foreach( test ${testSources} )

  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test} ${mock_object_headers} ${mock_object_impl}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${decoratedDependencies} gtest dataRepository )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )

endforeach()
