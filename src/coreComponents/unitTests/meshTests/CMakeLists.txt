# Specify list of tests
set( gtest_geosx_tests
     testMeshEnums.cpp
     testMeshGeneration.cpp
     testNeighborCommunicator.cpp )

set( gtest_geosx_mpi_tests
     testNeighborCommunicator.cpp )

if( ENABLE_VTK )
  list( APPEND gtest_geosx_tests
        testVTKImport.cpp )
  list( APPEND gtest_geosx_mpi_tests
        testVTKImport.cpp )
endif()

if( ENABLE_VTK )
  set( TEST_MESH_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/meshDirName.hpp.in ${CMAKE_BINARY_DIR}/include/tests/meshDirName.hpp )
endif()

set( dependencyList ${parallelDeps} gtest )

if( GEOSX_BUILD_SHARED_LIBS )
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

if (TARGET fmt::fmt-header-only)
  list( APPEND dependencyList fmt::fmt-header-only )
endif()

if (TARGET fmt)
  list( APPEND dependencyList fmt )
endif()

# Add gtest C++ based tests
foreach( test ${gtest_geosx_tests} )
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} ${CMAKE_CURRENT_LIST_DIR} )

  # Older g++ version need a hint to use std::filesystem
  if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
      if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9 )
          target_link_libraries( ${test_name} PRIVATE "stdc++fs" )
      endif ()
  endif ()
  if ( ENABLE_CUDA )
    target_link_libraries( ${test_name} PRIVATE "stdc++fs" )
  endif ()
endforeach()

if( ENABLE_MPI )

  set( nranks 2 )

  foreach( test ${gtest_geosx_mpi_tests} )
    get_filename_component( file_we ${test} NAME_WE )
    set( test_name ${file_we}_mpi )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
    # Older g++ version need a hint to use std::filesystem
    if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
        if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9 )
            target_link_libraries( ${test_name} PRIVATE "stdc++fs" )
        endif ()
    endif ()
    if ( ENABLE_CUDA )
      target_link_libraries( ${test_name} PRIVATE "stdc++fs" )
    endif ()
  endforeach()
endif()
