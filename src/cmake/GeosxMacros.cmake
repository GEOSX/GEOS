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
## combinatoricLists( RESULT <name of the variable to store the results in>
##                    JOIN <delimiter to use when combining items>
##                    LISTS <list of lists to generate combinations from>
##                    PREFIX <used for recursive expantsion, user should not supply this> )
##
## This function generates all combinations of elements from the provided lists.
## The results are combined into a single string, using the specified delimiter.
## The results are stored in the provided variable in the parent scope.
##
##------------------------------------------------------------------------------

function( combinatoricLists )
  set( options "" )
  set( oneValueArgs PREFIX RESULT JOIN )
  set( multiValueArgs LISTS )
  cmake_parse_arguments(PARSE_ARGV 0 ARG "${options}" "${oneValueArgs}" "${multiValueArgs}")

  # If PREFIX is not provided, set it to an empty string
  if(NOT DEFINED ARG_PREFIX)
    set(ARG_PREFIX "")
  endif()

  if(NOT DEFINED ARG_JOIN)
    set(ARG_JOIN "#")
  endif()

  if("${ARG_LISTS}" STREQUAL "")
    set(${ARG_RESULT} "${ARG_PREFIX}" PARENT_SCOPE)
  else()
    set(result "")
    list(GET ARG_LISTS 0 listname)
    list(REMOVE_AT ARG_LISTS 0)

    foreach(item IN LISTS ${listname})
      if( "${ARG_PREFIX}" STREQUAL "" )
        set( newprefix "${item}" )
      else()
        set( newprefix "${ARG_PREFIX}${ARG_JOIN}${item}" )
      endif()
      combinatoricLists( PREFIX "${newprefix}" LISTS ${ARG_LISTS} RESULT newresult )


      list( APPEND result ${newresult} )
    endforeach( )

    # Return results
    set( ${ARG_RESULT} ${result} PARENT_SCOPE )
  endif()
endfunction()

##------------------------------------------------------------------------------
## json_to_list( json_string <JSON array to convert>
##               result <name of the variable to store the results in> )
##
## This macro converts a JSON array into a CMake list.
## The result is stored in the provided variable.
##
##------------------------------------------------------------------------------
macro( json_to_list json_string result )
  # Get the size of the JSON array
  string ( JSON array_length LENGTH ${json_string} )

  # Initialize an empty list
  set ( cmake_list "" )

  if ( "${array_length}" GREATER 0 )

    # Decrement the array_length by 1 because the RANGE in CMake is inclusive
    math ( EXPR array_length "${array_length}-1" )

    # Loop over the JSON array and append each element to the CMake list
    foreach ( index RANGE ${array_length} )
      string ( JSON item GET ${json_string} ${index} )
      list ( APPEND cmake_list ${item} )
    endforeach ( )

  endif ( )

  # Set the result
  set ( ${result} ${cmake_list} )
endmacro ( )


##------------------------------------------------------------------------------
## json_list_of_pairs_to_two_lists( json_list <JSON list of pairs to convert>
##                                  first_elements_list <name of the variable to store the first elements in>
##                                  second_elements_list <name of the variable to store the second elements in> )
##
## This macro converts a JSON list of pairs into two CMake lists.
## The first elements of each pair are stored in the first result list,
## and the second elements are stored in the second result list.
##
##------------------------------------------------------------------------------
macro(json_list_of_pairs_to_two_lists json_list first_elements_list second_elements_list)
  # Initialize empty lists for the first and second elements
  set( first_elements "" )
  set( second_elements "" )

  # Get the length of the JSON list
  string ( JSON json_list_length LENGTH ${json_list} )

  if ( "${json_list_length}" GREATER 0 )
    # Decrement the array_length by 1 because the RANGE in CMake is inclusive
    math ( EXPR json_list_length "${json_list_length}-1" )

    # Iterate over the indices of the JSON list
    foreach( index RANGE 0 ${json_list_length} )
        # Get the sublist at the current index
        string( JSON sublist GET ${json_list} ${index} )

        # Get the elements of the sublist
        string( JSON first_element GET ${sublist} 0 )
        string( JSON second_element GET ${sublist} 1 )

        # Append the elements to their respective lists
        list( APPEND first_elements ${first_element} )
        list( APPEND second_elements ${second_element} )
    endforeach ( )
  endif ( )
  # Set the result in the parent scope
  set( ${first_elements_list} ${first_elements} )
  set( ${second_elements_list} ${second_elements} )
endmacro ( )

##------------------------------------------------------------------------------
## set_vars_from_lists( names <list of variable names>
##                      values <list of variable values> )
##
## This macro sets multiple variables at once.
## The names and values of the variables are provided as two lists.
## The variables are set in the current scope.
##
##------------------------------------------------------------------------------
macro(set_vars_from_lists names values)
  # Get the length of the lists
  list(LENGTH ${names} names_length)
  list(LENGTH ${values} values_length)

  # Check if the lengths of the lists are equal
  if(NOT ${names_length} EQUAL ${values_length})
    message(FATAL_ERROR "The lists 'names' and 'values' must be of the same length.")
  endif()

  # Iterate over the lists and set the variables
  math(EXPR length "${names_length} - 1")
  foreach(index RANGE ${length})
    list(GET ${names} ${index} name)
    list(GET ${values} ${index} value)
    set( ${name} ${value} )
  endforeach()
endmacro()

##------------------------------------------------------------------------------
## generateKernels( TEMPLATE <path to the template file>
##                  KEY <key to use when looking up the kernel specification in the JSON object>
##                  RESULT <name of the variable to store the list of generated files in>
##                  SPLIT <delimiter used to split the instantiation strings into lists>
##                  JSON <name of the JSON object variable> )
##
## This function generates kernel files from a template.
## The kernel specification is looked up in a JSON object using the provided key.
## The list of generated files is stored in the provided variable in the parent scope.
## The delimiter is used to split the instantiation strings into lists of symbol values.
## The JSON object is expected to be defined in the calling scope.
##
##------------------------------------------------------------------------------
function(generateKernels)
  set(options)
  set(oneValueArgs TEMPLATE KEY HEADERS SOURCES SPLIT JSON)
  set(multiValueArgs)

  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT DEFINED ARG_SPLIT)
    set(ARG_SPLIT "#")
  endif()

  # Get the relative directory of the file to the current source dir
  get_filename_component(localRelDir ${ARG_TEMPLATE} DIRECTORY)
  string(REPLACE ${CMAKE_SOURCE_DIR} "" relativeSourceDir ${CMAKE_CURRENT_SOURCE_DIR})

  string( JSON kernelSpec ERROR_VARIABLE jsonError GET ${${ARG_JSON}} ${ARG_KEY} )

  # Skip generation if key not found in kernel spec
  if(jsonError)
    message(FATAL_ERROR "${jsonError}, skipping kernel generation")
    return()
  endif()

  string( JSON symbolsJSON GET ${kernelSpec} "vars" )
  string( JSON constantsJSON GET ${kernelSpec} "constants" )
  string( JSON combinationsJSON GET ${kernelSpec} "combinations" )
  string( JSON explicitComboJSON GET ${kernelSpec} "explicit" )

  json_to_list( "${symbolsJSON}" symbolList )
  json_list_of_pairs_to_two_lists( "${constantsJSON}" constantSymbols constantValues )

  set( combinatoricSymbolList "" )
  string(LENGTH ${combinationsJSON} comboStrLen)
  if(comboStrLen GREATER 0)
    set(symbolListLists "")
    foreach(symbol IN LISTS symbolList)
      string(JSON symbolArray${symbol} GET ${combinationsJSON} ${symbol})
      json_to_list("${symbolArray${symbol}}" symbolList${symbol})
      list(APPEND symbolListLists symbolList${symbol})
    endforeach()
    combinatoricLists( LISTS ${symbolListLists}
                       JOIN ${ARG_SPLIT} # join with what we will later split using
                       RESULT combinatoricSymbolList )
  endif()
  json_to_list( "${explicitComboJSON}" explicitComboList )
  list( APPEND combinatoricSymbolList ${explicitComboList} )
  list( REMOVE_DUPLICATES combinatoricSymbolList )

  set(typeCombinationList "")
  set(generatedHeadersList "")
  set(generatedSourcesList "")

  list(LENGTH symbolList symbolCount)
  foreach(instantiation IN LISTS combinatoricSymbolList)
    string(REGEX REPLACE "[${ARG_SPLIT}\<\>\: ,.]" "_" suffix ${instantiation})
    set(generatedFileName "${CMAKE_BINARY_DIR}/generatedSrc/${relativeSourceDir}/${localRelDir}/${ARG_KEY}_${suffix}.cpp")
    string(REPLACE "${ARG_SPLIT}" ";" instantiationList ${instantiation})
    list(LENGTH instantiationList valueCount)
    if(symbolCount EQUAL valueCount)
      set_vars_from_lists( symbolList instantiationList )
      set_vars_from_lists( constantSymbols constantValues )
      message(STATUS "Generating file: ${generatedFileName}")
      configure_file(${ARG_TEMPLATE} ${generatedFileName} @ONLY)
      list(APPEND generatedSourcesList ${generatedFileName})

      string(REPLACE "${ARG_SPLIT}" ", " typeCombination ${instantiation})
      set(typeCombinationList "${typeCombinationList},
  types::TypeList< ${typeCombination} >")
    endif()
  endforeach()
  string(SUBSTRING "${typeCombinationList}" 4 -1 typeCombinationList)

  # Generate the dispatch type list header
  set(KERNEL_GROUP_NAME ${ARG_KEY})
  set(generatedFileName "${CMAKE_BINARY_DIR}/include/${relativeSourceDir}/${localRelDir}/${KERNEL_GROUP_NAME}DispatchTypeList.hpp")
  string(REPLACE "coreComponents/" "" generatedFileName "${generatedFileName}")
  message(STATUS "Generating file: ${generatedFileName}")
  configure_file(KernelDispatchTypeList.hpp.template "${generatedFileName}" @ONLY)
  list(APPEND generatedHeadersList ${generatedFileName})

  # The list of generated files is returned by setting the variable in the parent scope
  if(DEFINED ARG_HEADERS)
    set(${ARG_HEADERS} ${generatedHeadersList} PARENT_SCOPE)
  endif()
  if(DEFINED ARG_SOURCES)
    set(${ARG_SOURCES} ${generatedSourcesList} PARENT_SCOPE)
  endif()
endfunction()

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


macro( geos_decorate_link_dependencies )

    set( options )
    set( singleValueArgs LIST )
    set( multiValueArgs DEPENDENCIES )

    # Parse the arguments to the macro
    cmake_parse_arguments( arg
         "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    set( staticLibs "" )
    set( otherLibs "" )
    foreach( dep ${arg_DEPENDENCIES} )
        if( NOT TARGET ${dep} )
            message( FATAL_ERROR "Dependency ${dep} not found" )
        endif()

        get_target_property( targetType ${dep} TYPE)
        # message( "  ${dep} targetType = ${targetType}" )  # debug

        if (targetType STREQUAL STATIC_LIBRARY)
            list( APPEND staticLibs ${dep} )
        else()
            list( APPEND otherLibs ${dep} )
        endif()
    endforeach()

    # message( "  staticLibs = ${staticLibs}" )  # debug
    # message( "  otherLibs = ${otherLibs}" )  # debug

    string (REPLACE ";" "," staticLibsString "${staticLibs}")
    set( ${arg_LIST} "$<LINK_LIBRARY:WHOLE_ARCHIVE,${staticLibsString}>" ${otherLibs} )

endmacro( geos_decorate_link_dependencies )