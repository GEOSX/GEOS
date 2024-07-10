# Specify list of tests
set( gtest_geosx_tests
     testCapillaryPressure.cpp
     testCO2BrinePVTModels.cpp
     testCO2SpycherPruessModels.cpp
     testDamage.cpp
     testMultiFluidCO2Brine.cpp
     testMultiFluidDeadOil.cpp
     testMultiFluidLiveOil.cpp
     testRelPerm.cpp
     testRelPermHysteresis.cpp )

set( gtest_triaxial_xmls
     testTriaxial_druckerPragerExtended.xml
     testTriaxial_elasticIsotropic.xml
     testTriaxial_elasticIsotropicPressureDependent.xml
     testTriaxial_modifiedCamClay.xml
     testTriaxial_delftEggUseLinear.xml
     testTriaxial_modifiedCamClayVolumetric.xml
     testTriaxial_delftEggLoadPathDryUseLinear.xml
     testTriaxial_delftEggLoadPathWetUseLinear.xml
     testTriaxial_delftEggCase1.xml
     testTriaxial_delftEggCase2.xml )

set( gtest_pvt_xmls
     testPVT.xml
     testPVT_CO2Brine.xml
     testPVT_CO2BrineTables.xml
     testPVT_PhaseComposition.xml
   )

set( gtest_reactivefluid_xmls
     testReactiveFluid.xml )

set( dependencyList ${parallelDeps} gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  list( APPEND dependencyList geosx_core )
else()
  list( APPEND dependencyList ${geosx_core_libs} )
endif()

if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testMultiFluidCompositionalMultiphasePVTPackage.cpp )
    list( APPEND dependencyList PVTPackage )
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

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )

endforeach()

# Add triaxial xml based tests
blt_add_executable( NAME testTriaxial
                    SOURCES testTriaxial.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz" )

foreach(test ${gtest_triaxial_xmls})
    get_filename_component( test_name ${test} NAME_WE )
    geos_add_test( NAME ${test_name}
                   COMMAND testTriaxial -i ${CMAKE_CURRENT_LIST_DIR}/${test} )
endforeach()

# Add pvtdriver xml based tests
blt_add_executable( NAME testPVT
                    SOURCES testPVT.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz" )

foreach(test ${gtest_pvt_xmls})
    get_filename_component( test_name ${test} NAME_WE )
    geos_add_test( NAME ${test_name}
                   COMMAND testPVT -i ${CMAKE_CURRENT_LIST_DIR}/${test} )
endforeach()


# Add reactivefluiddriver xml based tests
blt_add_executable( NAME testReactiveFluid
                    SOURCES testReactiveFluid.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz" )

foreach(test ${gtest_reactivefluid_xmls})
    get_filename_component( test_name ${test} NAME_WE )
    geos_add_test( NAME ${test_name}
                   COMMAND testReactiveFluid -i ${CMAKE_CURRENT_LIST_DIR}/${test} )
endforeach()