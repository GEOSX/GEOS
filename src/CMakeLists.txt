
cmake_minimum_required( VERSION 3.23 )

# At the moment we are manually passing the cuda arch flag.
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.23.0")
  cmake_policy(SET CMP0104 OLD) # when using nvcc populate CMAKE_CUDA_ARCHITECTURES, raise error if we can't
endif()
cmake_policy(SET CMP0074 NEW) # dont ignore <PackageName>_ROOT env vars when searching for packages via find_package()
cmake_policy(SET CMP0066 NEW) # use CMAKE_<LANG>_FLAGS_<BUILD_TYPE> for try_compile() instead of only CMAKE_<LANG>_FLAGS
cmake_policy(SET CMP0056 NEW) # use CMAKE_EXE_LINKER_FLAGS in try_compile() in addition to CMAKE_<LANG>_FLAGS

################################
# GEOSX
################################
project( geosx LANGUAGES C CXX )
include(GNUInstallDirs)
set( BLT_CXX_STD "c++17" CACHE STRING "Version of C++ standard" )
if(CMAKE_BUILD_TYPE EQUAL "Debug")
  set( ENABLE_WARNINGS_AS_ERRORS "OFF" CACHE PATH "")
else()
  set( ENABLE_WARNINGS_AS_ERRORS "ON" CACHE PATH "")
endif()
set( CMAKE_ENABLE_EXPORTS ON )

# populated by `geos_add_test` (private list of all geos tests, keep in sync with the macro definition)
set_property( GLOBAL PROPERTY geos_tests_exe_list )

################################
# Include blt
################################
option( ENABLE_BENCHMARKS "Enables benchmarks" ON )
include( cmake/blt/SetupBLT.cmake )
set( BLT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/cmake/blt/ )

# use, i.e. don't skip the full RPATH for the build tree
set( CMAKE_SKIP_BUILD_RPATH FALSE )

# when building, don't use the install RPATH already
# (but later on when installing)
set( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
set( CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
endif()

# Removes -pthread from BLT CUDA libraries lists
if( ENABLE_CUDA )
  unset( _cuda_libraries )
  foreach( _lib ${BLT_CUDA_LIBRARIES} )
    if( "${_lib}" STREQUAL "-pthread" )
      list( APPEND _cuda_libraries "-lpthread" )
    else()
      list( APPEND _cuda_libraries "${_lib}" )
    endif()
  endforeach()

  set( BLT_CUDA_LIBRARIES "${_cuda_libraries}" CACHE STRING "" FORCE )
  set( BLT_CUDA_RUNTIME_LIBRARIES "${_cuda_libraries}" CACHE STRING "" FORCE )
endif()

################################
# Include standard build system logic and options
################################
include( cmake/CMakeBasics.cmake )
set( CMAKE_DEBUG_POSTFIX "" CACHE STRING "" FORCE )


################################
# Add Thirdparty Builtin Libs
################################
add_subdirectory( thirdparty )


################################
# Add components
################################
unset( externalComponentsLinkList )
add_subdirectory( coreComponents )
add_subdirectory( externalComponents )
add_subdirectory( main )


################################
# Add docs
################################
if ( ENABLE_DOCS )
  add_subdirectory( docs )
endif()


################################
# Deploy headers
################################
install( DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION . )


################################
# Generate version information
################################
include( cmake/GeosxVersion.cmake )
message( STATUS "Configuring GEOSX version ${GEOS_VERSION_FULL}" )


################################
# Create header of configuration options
################################
include( cmake/GeosxConfig.cmake )

################################
# Add main
################################
unset( extraComponentsLinkList )

if( ENABLE_OPENMP )
  list( APPEND extraComponentsLinkList openmp )
endif()

if( ENABLE_CUDA )
  list( APPEND extraComponentsLinkList cuda )
endif()

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND extraComponentsLinkList CUDA::nvToolsExt )
endif()

if( ENABLE_HIP )
  list( APPEND extraComponentsLinkList blt::hip )
endif()

if ( GEOS_BUILD_SHARED_LIBS )
  list( APPEND extraComponentsLinkList geosx_core )
else()
  list( APPEND extraComponentsLinkList ${geosx_core_libs} )
endif()

blt_add_executable( NAME geosx
                    SOURCES main/main.cpp
                    DEPENDS_ON ${extraComponentsLinkList}
                               ${externalComponentsLinkList} fmt::fmt-header-only HDF5::HDF5 )

# blt_print_target_properties(TARGET               geosx
#                             CHILDREN             TRUE
#                             PROPERTY_NAME_REGEX  "LINK_LIBRARIES"
#                             )

# blt_print_target_properties(TARGET               HDF5::HDF5
#                             CHILDREN             TRUE )

# Seems to be required on some CMake versions (e.g. 3.16) to get enforce device linking
if( ${ENABLE_HYPRE_DEVICE} STREQUAL "CUDA"  )
  set_target_properties( geosx PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS TRUE )
endif()


target_include_directories( geosx PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

# To change the runtime path during installation
set_target_properties( geosx PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
set_target_properties( geosx PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE )

install(TARGETS geosx RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
# Duplicate the installed `geosx` executable and rename it `geos` to let the end users adapt their scripts.
install( FILES ${CMAKE_INSTALL_PREFIX}/bin/geosx
         DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
         RENAME geos
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

if( ENABLE_XML_UPDATES AND ENABLE_MPI AND UNIX AND NOT CMAKE_HOST_APPLE AND NOT ENABLE_CUDA AND NOT ENABLE_HIP )

  set(SCHEMA_DIR ${CMAKE_SOURCE_DIR}/coreComponents/schema)
  set(SCRIPT_DIR ${CMAKE_SOURCE_DIR}/../scripts)

  add_custom_target( geosx_generate_schema
                     ALL
                     COMMAND geosx -s ${SCHEMA_DIR}/schema.xsd >generate_schema.log 2>&1 || (cat generate_schema.log && exit 1)
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx
                     COMMENT "Generating XML schema"
                     )

  add_custom_target( geosx_update_rst_tables
                     ALL
                     COMMAND python3 ${SCRIPT_DIR}/SchemaToRSTDocumentation.py ${SCHEMA_DIR} >update_rst_tables.log 2>&1 || (cat update_rst_tables.log && exit 1)
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx_generate_schema
                     COMMENT "Generating schema-dependent RST files"
                     )

  add_custom_target( geosx_validate_all_xml_files
                     ALL
                     COMMAND bash ${SCRIPT_DIR}/validateXMLFiles.bash -g ${SCHEMA_DIR}/schema.xsd ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../examples ${CMAKE_SOURCE_DIR}/../inputFiles
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx_update_rst_tables
                     COMMENT "Validating all XML files in the repository against the schema"
                     )
endif()

if( ENABLE_PYGEOSX )
   add_subdirectory( pygeosx )
endif()

if( SPHINX_FOUND )
  blt_add_sphinx_target( geosx_docs )
endif()


install( DIRECTORY ../examples/
         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${CMAKE_PROJECT_NAME}/examples )

install( DIRECTORY ../scripts/
         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${CMAKE_PROJECT_NAME}/scripts )

add_custom_target( geosx_generate_install_schema
         ALL
         COMMAND geosx -s ${CMAKE_BINARY_DIR}/schema.xsd >generate_schema.log 2>&1 || (cat generate_schema.log && exit 1)
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
         DEPENDS geosx
         COMMENT "Generating XML schema to install"
         )

if( NOT GEOS_INSTALL_SCHEMA )
  set_target_properties(geosx_generate_install_schema PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
endif()

install( FILES ${CMAKE_BINARY_DIR}/schema.xsd
         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${CMAKE_PROJECT_NAME}/schema OPTIONAL)


################################
# Add python environment setup
################################
# message(WARNING "Temporarily changing the geosPythonBranch to feature/packaging")
# set(GEOS_PYTHON_PACKAGES_BRANCH "feature/packaging" CACHE STRING "" FORCE)


if ( Python3_EXECUTABLE )
    message(STATUS "Found python version ${Python3_VERSION}")
    if (${Python3_VERSION} VERSION_LESS "3.6.0")
      message(STATUS "Note: try setting Python3_ROOT_DIR and/or Python3_EXECUTABLE in your host config to the appropriate version.")
      message( FATAL_ERROR "Building the GEOSX python tools requires Python >= 3.6." )
    endif()

    # Select the version of python to target
    if( ENABLE_PYGEOSX )
        set( PYTHON_POST_EXECUTABLE ${CMAKE_BINARY_DIR}/lib/PYGEOSX/bin/python CACHE PATH "" FORCE )
    else()
        set( PYTHON_POST_EXECUTABLE ${Python3_EXECUTABLE} CACHE PATH "" FORCE )
    endif()

    # Check for the virtualenv package
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import virtualenv"
        RESULT_VARIABLE VIRTUALENV_AVAILABLE
    )

    if (NOT ${VIRTUALENV_AVAILABLE} EQUAL 0)
        message(WARNING "The \"virtualenv\" package was not found in the target python environment (${Python3_EXECUTABLE}).  This package may be required to build PYGEOSX or the python development environment.")
    endif()

    # Check for the requested branch of geosPythonPackages (used for testing)
    if( NOT GEOS_PYTHON_PACKAGES_BRANCH )
      set(GEOS_PYTHON_PACKAGES_BRANCH "main" CACHE STRING "" FORCE)
    endif()

    set( GEOS_PYTHON_TOOLS_BINS
         "${CMAKE_BINARY_DIR}/bin/preprocess_xml"
         "${CMAKE_BINARY_DIR}/bin/format_xml" )

    add_custom_command( OUTPUT ${GEOS_PYTHON_TOOLS_BINS}
                        COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/setupPythonEnvironment.bash -p ${PYTHON_POST_EXECUTABLE} -b ${CMAKE_BINARY_DIR}/bin --python-pkg-branch ${GEOS_PYTHON_PACKAGES_BRANCH}
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )

    add_custom_target( geosx_python_tools
                       DEPENDS ${GEOS_PYTHON_TOOLS_BINS} )

    add_custom_target( geosx_python_tools_clean
                       COMMAND rm ${GEOS_PYTHON_TOOLS_BINS} )

    add_custom_target( geosx_python_tools_test
                       COMMAND ${CMAKE_BINARY_DIR}/python/geosx/bin/test_geosx_xml_tools
                       COMMAND rm -r ${CMAKE_BINARY_DIR}/python/geosx_xml_tools_tests*
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python
                       DEPENDS geosx_python_tools
                     )

    add_custom_target( geosx_format_all_xml_files
                       COMMAND bash ${CMAKE_SOURCE_DIR}/../scripts/formatXMLFiles.bash -g ${CMAKE_BINARY_DIR}/bin/format_xml ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../examples
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS geosx_xml_tools
                     )

else()
    message(WARNING "Building the GEOSX python tools requires Python >= 3.6.")
    message(STATUS "If you need these, try setting Python3_ROOT_DIR and/or Python3_EXECUTABLE in your host config.")
endif()


################################
# Add integratedTests
################################
if (NOT DEFINED ENABLE_ATS)
    set( ENABLE_ATS true CACHE BOOL "")
endif()


if ( ENABLE_ATS )
    if (NOT DEFINED ATS_WORKING_DIR)
        message( WARNING "ATS_WORKING_DIR is not defined (required for integrated testing system)" )
        message( WARNING "Defaulting to ${CMAKE_BINARY_DIR}/integratedTests/workingDir" )
        set( ATS_WORKING_DIR "${CMAKE_BINARY_DIR}/integratedTests/workingDir" CACHE PATH "")
    endif()

    if (NOT DEFINED ATS_BASELINE_DIR)
        message( WARNING "ATS_BASELINE_DIR is not defined (required for integrated testing system)" )
        message( WARNING "Defaulting to ${CMAKE_SOURCE_DIR}/../integratedTests" )
        set( ATS_BASELINE_DIR "${CMAKE_SOURCE_DIR}/../integratedTests" CACHE PATH "")
    endif()

    if (NOT Python3_EXECUTABLE)
        message( FATAL_ERROR "An appropriate version of python was not found (required for integrated testing system).  Try setting Python3_ROOT_DIR and/or Python3_EXECUTABLE in your host config." )
    endif() 

     # Setup testing
    set( ATS_SCRIPT
        "${CMAKE_BINARY_DIR}/integratedTests/geos_ats.sh"
       )

    add_custom_command( OUTPUT ${ATS_SCRIPT}
                        COMMAND ${CMAKE_BINARY_DIR}/bin/setup_ats_environment ${CMAKE_SOURCE_DIR}/.. ${CMAKE_BINARY_DIR} ${ATS_BASELINE_DIR} ${ATS_WORKING_DIR} ${ATS_ARGUMENTS}
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       )

    add_custom_target( ats_environment
                       DEPENDS geosx_python_tools
                       DEPENDS ${ATS_SCRIPT} )

    add_custom_target( ats_run
                       COMMAND ${CMAKE_BINARY_DIR}/integratedTests/geos_ats.sh
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS ats_environment
                     )

    add_custom_target( ats_clean
                       COMMAND ${CMAKE_BINARY_DIR}/integratedTests/geos_ats.sh -a veryclean
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS ats_environment
                     )

    add_custom_target( ats_rebaseline
                       COMMAND ${CMAKE_BINARY_DIR}/integratedTests/geos_ats.sh -a rebaseline
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS ats_environment
                     )

    add_custom_target( ats_rebaseline_failed
                       COMMAND ${CMAKE_BINARY_DIR}/integratedTests/geos_ats.sh -a rebaselinefailed
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       DEPENDS ats_environment
                     )
endif()


# Python formatting
if ( ENABLE_YAPF )
    set( integrated_tests_python_sources )
    file( GLOB_RECURSE integrated_tests_python_sources "${CMAKE_SOURCE_DIR}/../inputFiles/*.py" )
    set( integrated_tests_ats_sources )
    file( GLOB_RECURSE integrated_tests_ats_sources "${CMAKE_SOURCE_DIR}/../inputFiles/*.ats" )

    blt_add_code_checks( PREFIX  integrated_tests_yapf_style
                         SOURCES ${integrated_tests_python_sources} ${integrated_tests_ats_sources} ${CMAKE_SOURCE_DIR}/coreComponents/dummy.cpp
                         YAPF_CFG_FILE ${PROJECT_SOURCE_DIR}/yapf.cfg )
endif()


# the following adds a `build_test` CMake target such that running `$ make build_test test`
# builds the unit tests before running them
get_property( tmp GLOBAL PROPERTY geos_tests_exe_list )
add_custom_target( build_test COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${tmp} )

message( DEBUG "geos_tests: ${tmp} " )
