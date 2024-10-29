##------------------------------------------------------------------------------
## geosx_add_code_checks( PREFIX     <Prefix used for created targets>
##                        UNCRUSTIFY_CFG_FILE <path to style config file> 
##                        EXCLUDES   [path1 [path2 ...]])
##
## Adds code checks to all source files under this directory.
##
## PREFIX is used in the creation of all the underlying targets. For example:
## <PREFIX>_uncrustify_check.
##
## EXCLUDES is used to exclude any files from the code checks. It is done with
## a simple CMake reg exp MATCHES check.
##
##------------------------------------------------------------------------------
macro( geosx_add_code_checks )

    set( options )
    set( singleValueArgs PREFIX UNCRUSTIFY_CFG_FILE )
    set( multiValueArgs EXCLUDES )

    # Parse the arguments to the macro
    cmake_parse_arguments( arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    set(_all_sources )
    file( GLOB_RECURSE _all_sources
          "*.cpp" "*.hpp" "*.cxx" "*.hxx" "*.cc" "*.c" "*.h" "*.hh"
          "*.F" "*.f" "*.f90" "*.F90" )

    # Check for excludes
    if( NOT DEFINED arg_EXCLUDES )
        set( _sources ${_all_sources} )
    else()
        unset( _sources )
        foreach( _source ${_all_sources} )
            set( _to_be_excluded FALSE )
            foreach( _exclude ${arg_EXCLUDES} )
                if( ${_source} MATCHES ${_exclude} )
                    set( _to_be_excluded TRUE )
                    break()
                endif()
            endforeach()

            if( NOT ${_to_be_excluded} )
                list( APPEND _sources ${_source} )
            endif()
        endforeach()
    endif()

    if ( ENABLE_UNCRUSTIFY )
        blt_add_code_checks( PREFIX  ${arg_PREFIX}
                             SOURCES ${_sources}
                             UNCRUSTIFY_CFG_FILE ${PROJECT_SOURCE_DIR}/uncrustify.cfg )
    endif()

    if (ENABLE_COVERAGE)
        blt_add_code_coverage_target( NAME   ${arg_PREFIX}_coverage
                                      RUNNER ctest -E 'blt_gtest_smoke|blt_mpi_smoke|testUncrustifyCheck|testDoxygenCheck'
                                      SOURCE_DIRECTORIES ${PROJECT_SOURCE_DIR}/coreComponents )
    endif()

endmacro( geosx_add_code_checks )

##------------------------------------------------------------------------------
## geos_add_test( NAME       [name]
##                COMMAND    [command]
##                EXECUTABLE [executable] )
##
## Adds a test to the project, remaining arguments are forwarded to `blt_add_test`
## As a side effect, uses and populates a global property named `geos_tests_exe_list`
## initialized in `src/CMakeLists.txt`.
##------------------------------------------------------------------------------
macro( geos_add_test )

    set( options )
    set( singleValueArgs NAME EXECUTABLE )
    set( multiValueArgs COMMAND )

    # Parse the arguments to the macro
    cmake_parse_arguments( arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( NOT arg_EXECUTABLE )
      list( GET arg_COMMAND 0 _test_executable )
      set( arg_EXECUTABLE ${_test_executable} )
    endif()

    get_property( tmp GLOBAL PROPERTY geos_tests_exe_list )
    list( APPEND tmp ${arg_EXECUTABLE} )
    set_property( GLOBAL PROPERTY geos_tests_exe_list "${tmp}" )

    message( DEBUG "arg_NAME=${arg_NAME} arg_EXECUTABLE=${arg_EXECUTABLE} arg_COMMAND=${arg_COMMAND} ARGN=${ARGN}" )  # debug

    blt_add_test( NAME ${arg_NAME}
                  COMMAND ${arg_COMMAND} ${ARGN} )

endmacro( geos_add_test )
