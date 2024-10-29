# Specify list of tests
set( dependencyList gtest )

if ( GEOS_BUILD_SHARED_LIBS )
  list( APPEND dependencyList geosx_core )
else()
  list( APPEND dependencyList ${geosx_core_libs} )
endif()

set( newComponentTests
     testNewComponent.cpp )

# Add gtest C++ based tests
foreach(test ${newComponentTests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )
                        
    geos_add_test( NAME ${test_name} 
                   COMMAND ${test_name} )
endforeach()
