####################################
#
# 3rd Party Dependencies
#
# Setup all GEOS TPL
#
####################################


################################
# Helper macros & functions
################################

macro(find_and_import)
    set(singleValueArgs NAME HEADER)
    set(multiValueArgs INCLUDE_DIRECTORIES
                       LIBRARY_DIRECTORIES
                       LIBRARIES
                       EXTRA_LIBRARIES
                       DEPENDS )

    ## parse the arguments
    cmake_parse_arguments(arg
                          "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT DEFINED arg_NAME)
        message(FATAL_ERROR "The find_and_import required parameter NAME specifies the name of the library to import.")
    endif()

    if(NOT DEFINED arg_INCLUDE_DIRECTORIES)
        message(FATAL_ERROR "The find_and_import required parameter INCLUDE_DIRECTORIES specifies the directories to search for the given header.")
    endif()

    if(NOT DEFINED arg_LIBRARY_DIRECTORIES)
        message(FATAL_ERROR "The find_and_import required parameter LIBRARY_DIRECTORIES specifies the directories to search for the given libraries.")
    endif()

    if(NOT DEFINED arg_HEADER)
        message(FATAL_ERROR "The find_and_import required parameter HEADER specifies the header to search for.")
    endif()

    if(NOT DEFINED arg_LIBRARIES)
        message(FATAL_ERROR "The find_and_import required parameter LIBRARIES specifies the libraries to search for.")
    endif()

    find_path(${arg_NAME}_INCLUDE_DIR ${arg_HEADER}
              PATHS ${arg_INCLUDE_DIRECTORIES}
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    if(${arg_NAME}_INCLUDE_DIR STREQUAL ${arg_NAME}_INCLUDE_DIR-NOTFOUND)
        message(FATAL_ERROR "Could not find '${arg_HEADER}' in '${arg_INCLUDE_DIRECTORIES}'")
    endif()

    blt_find_libraries(FOUND_LIBS ${arg_NAME}_LIBRARIES
                       NAMES ${arg_LIBRARIES}
                       PATHS ${arg_LIBRARY_DIRECTORIES}
                       REQUIRED ON)

    blt_import_library(NAME ${arg_NAME}
                         INCLUDES ${${arg_NAME}_INCLUDE_DIR}
                         LIBRARIES ${${arg_NAME}_LIBRARIES} ${arg_EXTRA_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON
                         DEPENDS_ON ${arg_DEPENDS})

endmacro(find_and_import)


macro(extract_version_from_header)
    set(singleValueArgs NAME
                        HEADER
                        VERSION_STRING
                        MAJOR_VERSION_STRING
                        MINOR_VERSION_STRING
                        PATCH_VERSION_STRING )

    cmake_parse_arguments(arg
      "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    file(READ ${arg_HEADER} header_file)

    if(DEFINED arg_VERSION_STRING)

        if(DEFINED arg_VERSION_STRING AND "${header_file}" MATCHES "${arg_VERSION_STRING} *\"([^\"]*)\"")
            set(${arg_NAME}_VERSION "${CMAKE_MATCH_1}" CACHE STRING "" FORCE)
        endif()

    else()

        if(DEFINED arg_MAJOR_VERSION_STRING AND "${header_file}" MATCHES "${arg_MAJOR_VERSION_STRING} *([0-9]+)")
            set(ver_major "${CMAKE_MATCH_1}")
        endif()

        if(DEFINED arg_MINOR_VERSION_STRING AND "${header_file}" MATCHES "${arg_MINOR_VERSION_STRING} *([0-9]+)")
            set(ver_minor ".${CMAKE_MATCH_1}")
        endif()

        if(DEFINED arg_PATCH_VERSION_STRING AND "${header_file}" MATCHES "${arg_PATCH_VERSION_STRING} *([0-9]+)")
            set(ver_patch ".${CMAKE_MATCH_1}")
        endif()

        set(${arg_NAME}_VERSION "${ver_major}${ver_minor}${ver_patch}" CACHE STRING "" FORCE)

    endif()

    message(" ----> ${arg_NAME}_VERSION = ${${arg_NAME}_VERSION}")

endmacro( extract_version_from_header)


macro(mandatory_tpl_doesnt_exist
      CURRENT_TPL_NAME
      CURRENT_TPL_DIR_VAR)

    message(FATAL_ERROR
            "GEOSX requires ${CURRENT_TPL_NAME}, either :\n"
            "  - Verify that you provided a valid TPL installation directory (GEOSX_TPL_DIR = \"${GEOSX_TPL_DIR}\"),\n"
            "  - Or set ${CURRENT_TPL_DIR_VAR} to the ${CURRENT_TPL_NAME} installation directory (${CURRENT_TPL_DIR_VAR} = \"${${CURRENT_TPL_DIR_VAR}}\").\n")

endmacro(mandatory_tpl_doesnt_exist)


set(thirdPartyLibs "")


################################
# BLAS/LAPACK
################################
include(cmake/thirdparty/FindMathLibraries.cmake)

blt_import_library(NAME blas
                   TREAT_INCLUDES_AS_SYSTEM ON
                   LIBRARIES ${BLAS_LIBRARIES})

blt_import_library(NAME lapack
                   DEPENDS_ON blas
                   TREAT_INCLUDES_AS_SYSTEM ON
                   LIBRARIES ${LAPACK_LIBRARIES})

################################
# Intel MKL
################################
if(ENABLE_MKL)
    message(STATUS "Using Intel MKL")

    blt_import_library(NAME mkl
                         INCLUDES ${MKL_INCLUDE_DIRS}
                         LIBRARIES ${MKL_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(TRILINOS_DEPENDS mkl)
    set(thirdPartyLibs ${thirdPartyLibs} mkl)

################################
# IBM ESSL
################################
elseif(ENABLE_ESSL)
    message(STATUS "Using up IBM ESSL")

    blt_import_library(NAME essl
                         INCLUDES ${ESSL_INCLUDE_DIRS}
                         LIBRARIES ${ESSL_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(TRILINOS_DEPENDS essl)
    set(thirdPartyLibs ${thirdPartyLibs} essl)
else()
    set(TRILINOS_DEPENDS blas lapack)
    set(thirdPartyLibs ${thirdPartyLibs} blas lapack)
endif()

################################
# Conduit
################################
if(DEFINED CONDUIT_DIR)
    message(STATUS "CONDUIT_DIR = ${CONDUIT_DIR}")

    find_package(Conduit REQUIRED
                 PATHS ${CONDUIT_DIR}/lib/cmake
                 NO_DEFAULT_PATH)

    message( " ----> Conduit_VERSION = ${Conduit_VERSION}")


    set(CONDUIT_TARGETS conduit conduit_relay conduit_blueprint)
    foreach(targetName ${CONDUIT_TARGETS} )
        get_target_property(includeDirs
                            ${targetName}
                            INTERFACE_INCLUDE_DIRECTORIES)

        set_property(TARGET ${targetName}
                     APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                     ${includeDirs})
    endforeach()

    # Conduit uses our HDF5 and we need to propagate the above fix.
    # get_target_property(CONDUIT_RELAY_INTERFACE_INCLUDE_DIRECTORIES conduit_relay INTERFACE_INCLUDE_DIRECTORIES)
    # list(REMOVE_ITEM CONDUIT_RELAY_INTERFACE_INCLUDE_DIRECTORIES /usr/include)
    # set_target_properties(conduit_relay PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${CONDUIT_RELAY_INTERFACE_INCLUDE_DIRECTORIES})

    # get_target_property(CONDUIT_RELAY_INTERFACE_SYSTEM_INCLUDE_DIRECTORIES conduit_relay INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    # list(REMOVE_ITEM CONDUIT_RELAY_INTERFACE_SYSTEM_INCLUDE_DIRECTORIES /usr/include)
    # set_target_properties(conduit_relay PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${CONDUIT_RELAY_INTERFACE_SYSTEM_INCLUDE_DIRECTORIES})

    set(thirdPartyLibs ${thirdPartyLibs} conduit::conduit)
else()
    mandatory_tpl_doesnt_exist("Conduit" CONDUIT_DIR)
endif()

################################
# HDF5
################################
if(DEFINED HDF5_DIR)
    message(STATUS "HDF5_DIR = ${HDF5_DIR}")

    set(HDF5_ROOT ${HDF5_DIR})
    set(HDF5_USE_STATIC_LIBRARIES FALSE)
    set(HDF5_NO_FIND_PACKAGE_CONFIG_FILE ON)
    include(FindHDF5)

    # On some platforms (Summit) HDF5 lists /usr/include in it's list of include directories.
    # When this happens you can get really opaque include errors.
    list(REMOVE_ITEM HDF5_INCLUDE_DIRS /usr/include)

    blt_import_library(NAME hdf5
                       INCLUDES ${HDF5_INCLUDE_DIRS}
                       LIBRARIES ${HDF5_LIBRARIES}
                       TREAT_INCLUDES_AS_SYSTEM ON)

    file(READ "${HDF5_DIR}/include/H5public.h" header_file )
    string(REGEX MATCH "version: *([0-9]+.[0-9]+.[0-9]+)" _ ${header_file})
    set( HDF5_VERSION "${CMAKE_MATCH_1}" CACHE STRING "" FORCE )
    message( " ----> HDF5 version ${HDF5_VERSION}")

    set(ENABLE_HDF5 ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} hdf5)
else()
    mandatory_tpl_doesnt_exist("hdf5" HDF5_DIR)
endif()

################################
# SILO
################################
if(DEFINED SILO_DIR AND ENABLE_SILO)
    message(STATUS "SILO_DIR = ${SILO_DIR}")

    find_and_import(NAME silo
                      INCLUDE_DIRECTORIES ${SILO_DIR}/include
                      LIBRARY_DIRECTORIES ${SILO_DIR}/lib
                      HEADER silo.h
                      LIBRARIES siloh5
                      DEPENDS hdf5)


    set(ENABLE_SILO ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} silo)
else()
    message(STATUS "Not using SILO.")
endif()

################################
# PUGIXML
################################
if(DEFINED PUGIXML_DIR)
    message(STATUS "PUGIXML_DIR = ${PUGIXML_DIR}")

    set(ENABLE_PUGIXML ON CACHE BOOL "")

    find_package(pugixml REQUIRED
                 PATHS ${PUGIXML_DIR}
                 NO_DEFAULT_PATH)

    message( " ----> pugixml_VERSION = ${pugixml_VERSION}")

    if(TARGET pugixml::pugixml)
      set(thirdPartyLibs ${thirdPartyLibs} pugixml::pugixml)
    endif()
    if(TARGET pugixml)
      set(thirdPartyLibs ${thirdPartyLibs} pugixml)
    endif()
else()
    mandatory_tpl_doesnt_exist("pugixml" PUGIXML_DIR)
endif()

################################
# CUDA
################################
if ( ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
  message( " ----> $CUDAToolkit_VERSION = ${CUDAToolkit_VERSION}")
endif()

################################
# CAMP ( required before raja on crusher / using spack installed tpls )
################################
if(DEFINED CAMP_DIR)
    if( CAMP_STANDALONE )
        # Should be found by raja, but it is possible for spack to misconfig raja so we need to find it
        message(STATUS "CAMP_DIR = ${CAMP_DIR}")
        find_package(camp REQUIRED PATHS ${CAMP_DIR} NO_DEFAULT_PATH)
        get_target_property(CAMP_INCLUDE_DIRS camp INTERFACE_INCLUDE_DIRECTORIES)
        set_target_properties(camp PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CAMP_INCLUDE_DIRS}")
    endif( )
endif()

################################
# RAJA
################################
if(DEFINED RAJA_DIR)
    message(STATUS "RAJA_DIR = ${RAJA_DIR}")
    find_package(RAJA REQUIRED
                 PATHS ${RAJA_DIR}
                 NO_DEFAULT_PATH)

    message( " ----> RAJA_VERSION = ${RAJA_VERSION}")

    get_target_property(RAJA_INCLUDE_DIRS RAJA INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(RAJA PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${RAJA_INCLUDE_DIRS}")
    set(ENABLE_RAJA ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} RAJA )
else()
    mandatory_tpl_doesnt_exist("RAJA" RAJA_DIR)
endif()

################################
# CAMP ( required after raja on lassen / using non-spack installed tpls )
################################
if(DEFINED CAMP_DIR)
    if( NOT DEFINED CAMP_STANDALONE OR NOT CAMP_STANDALONE )
        # Should be found by raja, but it is possible for spack to misconfig raja so we need to find it
        message(STATUS "CAMP_DIR = ${CAMP_DIR}")
        find_package(camp REQUIRED PATHS ${CAMP_DIR} NO_DEFAULT_PATH)
        get_target_property(CAMP_INCLUDE_DIRS camp INTERFACE_INCLUDE_DIRECTORIES)
        set_target_properties(camp PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CAMP_INCLUDE_DIRS}")
    endif()
endif()

################################
# Umpire
################################
if(DEFINED UMPIRE_DIR)
    message(STATUS "UMPIRE_DIR = ${UMPIRE_DIR}")

    find_package(umpire REQUIRED
                 PATHS ${UMPIRE_DIR}
                 NO_DEFAULT_PATH)

    message( " ----> umpire_VERSION = ${umpire_VERSION}")

    set(ENABLE_UMPIRE ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} umpire)
else()
    mandatory_tpl_doesnt_exist("Umpire" UMPIRE_DIR)
endif()


################################
# CHAI
################################
if(DEFINED CHAI_DIR)
    message(STATUS "CHAI_DIR = ${CHAI_DIR}")

    find_package(chai REQUIRED
                 PATHS ${CHAI_DIR}
                 NO_DEFAULT_PATH)

    message( " ----> chai_VERSION = ${chai_VERSION}")

    get_target_property(CHAI_INCLUDE_DIRS chai INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(chai
                          PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CHAI_INCLUDE_DIRS}")

    set(ENABLE_CHAI ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} chai)
else()
    mandatory_tpl_doesnt_exist("CHAI" CHAI_DIR)
endif()

################################
# Adiak
################################
if(DEFINED ADIAK_DIR)
    message(STATUS "ADIAK_DIR = ${ADIAK_DIR}")

    find_package(adiak REQUIRED
                 PATHS ${ADIAK_DIR}
                 NO_DEFAULT_PATH)

    # Header file provides incorrect version 0.3.0
    message( " ----> adiak_VERSION = 0.2.2")

    set(adiak_target "")
    if(TARGET adiak::adiak)
      set(adiak_target ${adiak_target} adiak::adiak)
    endif()
    if(TARGET adiak)
      set(adiak_target ${adiak_target} adiak)
    endif()

    set_property(TARGET ${adiak_target}
                 APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                 ${adiak_INCLUDE_DIR} )
    set_property(TARGET ${adiak_target}
                 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                 ${adiak_INCLUDE_DIR} )

    set(ENABLE_ADIAK ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} ${adiak_target})
else()
    if(ENABLE_ADIAK)
        message(WARNING "ENABLE_ADIAK is ON but ADIAK_DIR isn't defined.")
    endif()

    set(ENABLE_ADIAK OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Adiak.")
endif()

################################
# Caliper
################################
if(DEFINED CALIPER_DIR)
    message(STATUS "CALIPER_DIR = ${CALIPER_DIR}")

    find_package(caliper REQUIRED
                 PATHS ${CALIPER_DIR}
                 NO_DEFAULT_PATH)

    extract_version_from_header( NAME caliper
                                 HEADER "${CALIPER_DIR}/include/caliper/caliper-config.h"
                                 MAJOR_VERSION_STRING "CALIPER_MAJOR_VERSION"
                                 MINOR_VERSION_STRING "CALIPER_MINOR_VERSION"
                                 PATCH_VERSION_STRING "CALIPER_PATCH_VERSION")

    set_property(TARGET caliper
                 APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                 ${caliper_INCLUDE_PATH} )

    set_property(TARGET caliper
                 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                 ${caliper_INCLUDE_PATH} )

    set(ENABLE_CALIPER ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} caliper)
else()
    if(ENABLE_CALIPER)
        message(WARNING "ENABLE_CALIPER is ON but CALIPER_DIR isn't defined.")
    endif()

    set(ENABLE_CALIPER OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Caliper.")
endif()

################################
# MATHPRESSO
################################
if(DEFINED MATHPRESSO_DIR)
    message(STATUS "MATHPRESSO_DIR = ${MATHPRESSO_DIR}")

    find_and_import(NAME mathpresso
                      INCLUDE_DIRECTORIES ${MATHPRESSO_DIR}/include
                      LIBRARY_DIRECTORIES ${MATHPRESSO_DIR}/lib
                      HEADER mathpresso/mathpresso.h
                      LIBRARIES mathpresso)

    set(ENABLE_MATHPRESSO ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} mathpresso)
else()
    if(ENABLE_MATHPRESSO)
        message(WARNING "ENABLE_MATHPRESSO is ON but MATHPRESSO_DIR isn't defined.")
    endif()

    set(ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using mathpresso.")
endif()

################################
# METIS
################################
if(DEFINED METIS_DIR)
    message(STATUS "METIS_DIR = ${METIS_DIR}")

    find_and_import(NAME metis
                      INCLUDE_DIRECTORIES ${METIS_DIR}/include
                      LIBRARY_DIRECTORIES ${METIS_DIR}/lib
                      HEADER metis.h
                      LIBRARIES metis)

    extract_version_from_header( NAME metis
                                 HEADER "${METIS_DIR}/include/metis.h"
                                 MAJOR_VERSION_STRING "METIS_VER_MAJOR"
                                 MINOR_VERSION_STRING "METIS_VER_MINOR"
                                 PATCH_VERSION_STRING "METIS_VER_PATCH")


    set(ENABLE_METIS ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} metis)
else()
    if(ENABLE_METIS)
        message(WARNING "ENABLE_METIS is ON but METIS_DIR isn't defined.")
    endif()

    set(ENABLE_METIS OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using METIS.")
endif()

################################
# PARMETIS
################################
if(DEFINED PARMETIS_DIR)
    message(STATUS "PARMETIS_DIR = ${PARMETIS_DIR}")

    find_and_import(NAME parmetis
                      INCLUDE_DIRECTORIES ${PARMETIS_DIR}/include
                      LIBRARY_DIRECTORIES ${PARMETIS_DIR}/lib
                      HEADER parmetis.h
                      LIBRARIES parmetis
                      DEPENDS metis)

    extract_version_from_header( NAME parmetis
                                 HEADER "${PARMETIS_DIR}/include/parmetis.h"
                                 MAJOR_VERSION_STRING "PARMETIS_MAJOR_VERSION"
                                 MINOR_VERSION_STRING "PARMETIS_MINOR_VERSION"
                                 PATCH_VERSION_STRING "PARMETIS_PATCH_VERSION")

    set(ENABLE_PARMETIS ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} parmetis)
else()
    if(ENABLE_PARMETIS)
        message(WARNING "ENABLE_PARMETIS is ON but PARMETIS_DIR isn't defined.")
    endif()

    set(ENABLE_PARMETIS OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using ParMETIS.")
endif()

################################
# SCOTCH
################################
if(DEFINED SCOTCH_DIR)
    message(STATUS "SCOTCH_DIR = ${SCOTCH_DIR}")

    find_and_import(NAME scotch
                      INCLUDE_DIRECTORIES ${SCOTCH_DIR}/include
                      LIBRARY_DIRECTORIES ${SCOTCH_DIR}/lib
                      HEADER scotch.h
                      LIBRARIES scotch scotcherr )

    find_and_import(NAME ptscotch
                      INCLUDE_DIRECTORIES ${SCOTCH_DIR}/include
                      LIBRARY_DIRECTORIES ${SCOTCH_DIR}/lib
                      DEPENDS scotch
                      HEADER ptscotch.h
                      LIBRARIES ptscotch ptscotcherr )

    extract_version_from_header( NAME scotch
                                 HEADER "${SCOTCH_DIR}/include/scotch.h"
                                 MAJOR_VERSION_STRING "SCOTCH_VERSION"
                                 MINOR_VERSION_STRING "SCOTCH_RELEASE"
                                 PATCH_VERSION_STRING "SCOTCH_PATCHLEVEL")

    set(ENABLE_SCOTCH ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} scotch ptscotch)
else()
    if(ENABLE_SCOTCH)
        message(WARNING "ENABLE_SCOTCH is ON but SCOTCH_DIR isn't defined.")
    endif()

    set(ENABLE_SCOTCH OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using SCOTCH.")
endif()

################################
# SUPERLU_DIST
################################
if(DEFINED SUPERLU_DIST_DIR)
    message(STATUS "SUPERLU_DIST_DIR = ${SUPERLU_DIST_DIR}")

    find_and_import(NAME superlu_dist
                      INCLUDE_DIRECTORIES ${SUPERLU_DIST_DIR}/include
                      LIBRARY_DIRECTORIES ${SUPERLU_DIST_DIR}/lib PATHS ${SUPERLU_DIST_DIR}/lib64
                      HEADER superlu_defs.h
                      LIBRARIES superlu_dist
                      DEPENDS parmetis blas lapack)


    extract_version_from_header( NAME superlu_dist
                                 HEADER "${SUPERLU_DIST_DIR}/include/superlu_defs.h"
                                 MAJOR_VERSION_STRING "SUPERLU_DIST_MAJOR_VERSION"
                                 MINOR_VERSION_STRING "SUPERLU_DIST_MINOR_VERSION"
                                 PATCH_VERSION_STRING "SUPERLU_DIST_PATCH_VERSION")

    set(ENABLE_SUPERLU_DIST ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} superlu_dist)
else()
    if(ENABLE_SUPERLU_DIST)
        message(WARNING "ENABLE_SUPERLU_DIST is ON but SUPERLU_DIST_DIR isn't defined.")
    endif()

    set(ENABLE_SUPERLU_DIST OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using superlu_dist.")
endif()

################################
# SUITESPARSE
################################
if(DEFINED SUITESPARSE_DIR)
    message(STATUS "SUITESPARSE_DIR = ${SUITESPARSE_DIR}")

    find_and_import(NAME suitesparse
                      INCLUDE_DIRECTORIES ${SUITESPARSE_DIR}/include
                      LIBRARY_DIRECTORIES ${SUITESPARSE_DIR}/lib ${SUITESPARSE_DIR}/lib64
                      HEADER umfpack.h
                      LIBRARIES umfpack
                      DEPENDS blas lapack)

    extract_version_from_header( NAME suitesparse
                                 HEADER "${SUITESPARSE_DIR}/include/umfpack.h"
                                 MAJOR_VERSION_STRING "UMFPACK_MAIN_VERSION"
                                 MINOR_VERSION_STRING "UMFPACK_SUB_VERSION"
                                 PATCH_VERSION_STRING "UMFPACK_SUBSUB_VERSION")

    set(ENABLE_SUITESPARSE ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} suitesparse)
else()
    if(ENABLE_SUITESPARSE)
        message(WARNING "ENABLE_SUITESPARSE is ON but SUITESPARSE_DIR isn't defined.")
    endif()

    set(ENABLE_SUITESPARSE OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using SuiteSparse.")
endif()

################################
# HYPRE
################################
if(DEFINED HYPRE_DIR AND ENABLE_HYPRE)
    message(STATUS "HYPRE_DIR = ${HYPRE_DIR}")

    set( HYPRE_DEPENDS blas lapack umpire )
    if( ENABLE_SUPERLU_DIST )
        list( APPEND HYPRE_DEPENDS superlu_dist )
    endif()
    if( ${ENABLE_HYPRE_DEVICE} STREQUAL "CUDA" )
        list( APPEND HYPRE_DEPENDS CUDA::cusparse CUDA::cublas CUDA::curand CUDA::cusolver )

        # Add libnvJitLink when using CUDA >= 12.2.2. Note: requires cmake >= 3.26
        if( CUDAToolkit_VERSION VERSION_GREATER_EQUAL "12.2.2" )
           list( APPEND HYPRE_DEPENDS CUDA::nvJitLink )
        endif()
    elseif( ${ENABLE_HYPRE_DEVICE} STREQUAL "HIP" )
        find_package( rocblas REQUIRED )
        find_package( rocsolver REQUIRED )
        find_package( rocsparse REQUIRED )
        find_package( rocrand REQUIRED )
        append( APPEND HYPRE_DEPENDS roc::rocblas roc::rocsparse roc::rocsolver roc::rocrand )
    endif( )

    find_and_import( NAME hypre
                     INCLUDE_DIRECTORIES ${HYPRE_DIR}/include
                     LIBRARY_DIRECTORIES ${HYPRE_DIR}/lib
                     HEADER HYPRE.h
                     LIBRARIES HYPRE
                     DEPENDS ${HYPRE_DEPENDS} )

    extract_version_from_header( NAME hypre
                                 HEADER "${HYPRE_DIR}/include/HYPRE_config.h"
                                 VERSION_STRING "HYPRE_RELEASE_VERSION" )

    # Extract some additional information about development version of hypre
    file( READ ${HYPRE_DIR}/include/HYPRE_config.h header_file )
    if( "${header_file}" MATCHES "HYPRE_DEVELOP_STRING *\"([^\"]*)\"" )
        set( hypre_dev_string "${CMAKE_MATCH_1}" )
        if( "${header_file}" MATCHES "HYPRE_BRANCH_NAME *\"([^\"]*)\"" )
            set( hypre_dev_branch "${CMAKE_MATCH_1}" )
        endif()
        set( hypre_VERSION "${hypre_dev_string} (${hypre_dev_branch})" CACHE STRING "" FORCE )
        message( " ----> hypre_VERSION = ${hypre_VERSION}" )
    endif()

    # Prepend Hypre to link flags, fix for Umpire appearing before Hypre on the link line
    # if (NOT CMAKE_HOST_APPLE)
    #   blt_add_target_link_flags (TO hypre FLAGS "-Wl,--whole-archive ${HYPRE_DIR}/lib/libHYPRE.a -Wl,--no-whole-archive")
    # endif()

    # if( ENABLE_CUDA AND ( NOT ${ENABLE_HYPRE_DEVICE} STREQUAL "CUDA" ) )
    #   set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)
    #   if( GEOSX_LA_INTERFACE STREQUAL "Hypre")
    #     message( FATAL_ERROR "Hypre LAI selected, but ENABLE_HYPRE_DEVICE not 'CUDA' while ENABLE_CUDA is ON.")
    #   endif()
    # else()
    #   set(ENABLE_HYPRE ON CACHE BOOL "")
    # endif()

    set( ENABLE_HYPRE ON CACHE BOOL "" )
    set( thirdPartyLibs ${thirdPartyLibs} hypre ${HYPRE_DEPENDS} )
else()
    if(ENABLE_HYPRE)
        message(WARNING "ENABLE_HYPRE is ON but HYPRE_DIR isn't defined.")
    endif()

    set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using HYPRE.")
endif()

################################
# TRILINOS
################################
if(DEFINED TRILINOS_DIR AND ENABLE_TRILINOS)
    message(STATUS "TRILINOS_DIR = ${TRILINOS_DIR}")

    include(${TRILINOS_DIR}/lib/cmake/Trilinos/TrilinosConfig.cmake)

    list(REMOVE_ITEM Trilinos_LIBRARIES "gtest")
    list(REMOVE_DUPLICATES Trilinos_LIBRARIES)

    blt_import_library(NAME trilinos
                       DEPENDS_ON ${TRILINOS_DEPENDS}
                       INCLUDES ${Trilinos_INCLUDE_DIRS}
                       LIBRARIES ${Trilinos_LIBRARIES}
                       TREAT_INCLUDES_AS_SYSTEM ON)

    extract_version_from_header( NAME trilinos
                                 HEADER "${TRILINOS_DIR}/include/Trilinos_version.h"
                                 VERSION_STRING "TRILINOS_VERSION_STRING" )

    # This conditional is due to the lack of mixedInt support on hypre GPU.
    # This can be removed when support is added into hypre.
    if( NOT ${ENABLE_HYPRE_DEVICE} STREQUAL "HIP" )
        set(ENABLE_TRILINOS ON CACHE BOOL "")
    endif()
    set(thirdPartyLibs ${thirdPartyLibs} trilinos)
else()
    if(ENABLE_TRILINOS)
        message(WARNING "ENABLE_TRILINOS is ON but TRILINOS_DIR isn't defined.")
    endif()

    set(ENABLE_TRILINOS OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Trilinos")
endif()

###############################
# PETSC
###############################
if(DEFINED PETSC_DIR AND ENABLE_PETSC)
    message(STATUS "PETSC_DIR = ${PETSC_DIR}")

    set( PETSC_DEPENDS metis blas lapack )
    if( ${ENABLE_SUPERLU_DIST} )
        set( PETSC_DEPENDS ${PETSC_DEPENDS} superlu_dist )
    endif()

    find_and_import(NAME petsc
                      INCLUDE_DIRECTORIES ${PETSC_DIR}/include
                      LIBRARY_DIRECTORIES ${PETSC_DIR}/lib
                      HEADER petscvec.h
                      LIBRARIES petsc
                      DEPENDS ${PETSC_DEPENDS})

    extract_version_from_header( NAME petsc
                                 HEADER "${PETSC_DIR}/include/petscversion.h"
                                 MAJOR_VERSION_STRING "PETSC_VERSION_MAJOR"
                                 MINOR_VERSION_STRING "PETSC_VERSION_MINOR"
                                 PATCH_VERSION_STRING "PETSC_VERSION_SUBMINOR")

    set(ENABLE_PETSC ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} petsc)
else()
    if(ENABLE_PETSC)
        message(WARNING "ENABLE_PETSC is ON but PETSC_DIR isn't defined.")
    endif()

    set(ENABLE_PETSC OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using PETSc")
endif()
################################
# VTK
################################
if(DEFINED VTK_DIR)
    message(STATUS "VTK_DIR = ${VTK_DIR}")
    find_package(VTK REQUIRED
                 PATHS ${VTK_DIR}
                 NO_DEFAULT_PATH)

    message( " ----> VTK_VERSION=${VTK_VERSION}")

    set( VTK_TARGETS
         VTK::FiltersParallelDIY2
         VTK::IOLegacy
         VTK::IOParallelXML
         VTK::IOXML
         VTK::ParallelMPI
         )
    foreach( targetName ${VTK_TARGETS} )

        get_target_property( includeDirs ${targetName}  INTERFACE_INCLUDE_DIRECTORIES)

        set_property(TARGET ${targetName}
                     APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                     ${includeDirs} )
    endforeach()

    set(ENABLE_VTK ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} VTK)
else()
    if(ENABLE_VTK)
        message(WARNING "ENABLE_VTK is ON but VTK_DIR isn't defined.")
    endif()

    set(ENABLE_VTK OFF CACHE BOOL "")
    message(STATUS "Not using VTK")
endif()

################################
# FMT
################################
if(DEFINED FMT_DIR)
    message(STATUS "FMT_DIR = ${FMT_DIR}")

    find_package(fmt REQUIRED
                 PATHS ${FMT_DIR}
                 NO_DEFAULT_PATH)

    message( " ----> fmt_VERSION = ${fmt_VERSION}")

    get_target_property(includeDirs fmt::fmt INTERFACE_INCLUDE_DIRECTORIES)

    set_property(TARGET fmt::fmt
                 APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                 ${includeDirs})

    set(ENABLE_FMT ON CACHE BOOL "")

    set(thirdPartyLibs ${thirdPartyLibs} fmt::fmt )
else()
    mandatory_tpl_doesnt_exist("{fmt}" FMT_DIR)
endif()

################################
# uncrustify
################################
if(UNCRUSTIFY_FOUND)
    message(STATUS "UNCRUSTIFY_EXECUTABLE = ${UNCRUSTIFY_EXECUTABLE}")

    # targets for verifying formatting
    if( NOT TARGET uncrustify_check )
        add_custom_target(uncrustify_check)
        add_dependencies(check uncrustify_check)
    endif()

    # targets for modifying formatting
    if( NOT TARGET uncrustify_style )
        add_custom_target(uncrustify_style)
        add_dependencies(style uncrustify_style)
    endif()
else()
    message(STATUS "Not using uncrustify.")
endif()

################################
# doxygen
################################
if(DOXYGEN_FOUND)
    message(STATUS "DOXYGEN_EXECUTABLE = ${DOXYGEN_EXECUTABLE}")
else()
    message(STATUS "Not using doxygen.")
endif()

################################
# Python
################################
message(" CMAKE_VERSION ${CMAKE_VERSION} ")
if( ${CMAKE_VERSION} VERSION_LESS "3.19" )
    set( PYTHON_AND_VERSION Python3 )
    set( PYTHON_OPTIONAL_COMPONENTS)
else()
    set( PYTHON_AND_VERSION Python3 3.6.0...3.12.2 )
    set( PYTHON_OPTIONAL_COMPONENTS OPTIONAL_COMPONENTS Development NumPy)
endif()
if(ENABLE_PYGEOSX)
    find_package(${PYTHON_AND_VERSION} REQUIRED
                 COMPONENTS Development NumPy)

    message( " ----> $Python3_VERSION = ${Python3_VERSION}")

    message(STATUS "Python3_EXECUTABLE=${Python3_EXECUTABLE}")
    message(STATUS "Python3_INCLUDE_DIRS = ${Python3_INCLUDE_DIRS}")
    message(STATUS "Python3_LIBRARY_DIRS = ${Python3_LIBRARY_DIRS}")
    message(STATUS "Python3_NumPy_INCLUDE_DIRS = ${Python3_NumPy_INCLUDE_DIRS}")

    if(DEFINED ENABLE_PYLVARRAY AND NOT ENABLE_PYLVARRAY)
        message(FATAL_ERROR "Cannot build pygeosx without pylvarray")
    else()
        set(ENABLE_PYLVARRAY ON CACHE BOOL "")
    endif()

    set(thirdPartyLibs ${thirdPartyLibs} Python3::Python Python3::NumPy)
else()
    message(STATUS "Not building pygeosx.")
    find_package(${PYTHON_AND_VERSION} ${PYTHON_OPTIONAL_COMPONENTS})
    message(STATUS "Python3_EXECUTABLE=${Python3_EXECUTABLE}")
endif()

################################
# LAI
################################
string(TOUPPER "${GEOSX_LA_INTERFACE}" upper_LAI)
if(NOT ENABLE_${upper_LAI})
  message(FATAL_ERROR "${GEOSX_LA_INTERFACE} LA interface is selected, but ENABLE_${upper_LAI} is OFF")
endif()
option(GEOSX_LA_INTERFACE_${upper_LAI} "${upper_LAI} LA interface is selected" ON)


message(STATUS "thirdPartyLibs = ${thirdPartyLibs}")

###############################
# NvToolExt
###############################
if ( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  set(thirdPartyLibs ${thirdPartyLibs} CUDA::nvToolsExt)
endif()

message(STATUS "thirdPartyLibs = ${thirdPartyLibs}")
