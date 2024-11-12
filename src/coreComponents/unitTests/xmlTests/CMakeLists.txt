# Specify list of tests
set( gtest_geosx_tests
     testXML.cpp )

set( gtest_tests_with_input
     testXMLFile.cpp )

set( input_files
     basic_input
     include_input )

set( multi_files
     multifile_input/fieldspec.xml
     multifile_input/mesh.xml
     multifile_input/outputs.xml
     multifile_input/solver.xml )

set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

# create test executables
foreach(test ${gtest_geosx_tests} ${gtest_tests_with_input})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} "-lz" )
endforeach()

# add tests that don't require input files
foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()

# add tests that do require input files
foreach( test ${gtest_tests_with_input} )
  get_filename_component( test_name ${test} NAME_WE )

  # test a single input file
  foreach( input ${input_files} )
    geos_add_test( NAME ${test_name}_${input}
                   EXECUTABLE ${test_name}
                   COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/${input}.xml )
  endforeach()

  # test multiple input files
  foreach( file ${multi_files} )
    set( args ${args} -i ${CMAKE_CURRENT_LIST_DIR}/${file} )
  endforeach()
  geos_add_test( NAME ${test_name}_multiple_input
                 COMMAND ${test_name} ${args} )
endforeach()
