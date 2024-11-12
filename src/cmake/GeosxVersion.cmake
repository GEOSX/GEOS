# This script generates version info for GEOSX from VERSION file and git.
# At config time, it sets various version variables and generates a build
# target named 'generate_version', that the code must set a dependency on.
# The target invokes this script again at build time, where it generates a
# version header with the most up-to-date version info. Thus any change to
# VERSION file or git state (in particular, branch name or HEAD commit sha)
# will trigger a partial rebuild (anything that includes version header).
#
# This self-referential approach inspired in part by:
# https://github.com/andrew-hardin/cmake-git-version-tracking

# In script mode CMAKE_SOURCE_DIR is set to the build directory,
# therefore the script uses SOURCE_DIR as the CMake root path.
if( NOT DEFINED CMAKE_SCRIPT_MODE_FILE )
  set( SOURCE_DIR ${CMAKE_SOURCE_DIR} )
endif()

# Get GEOSX version from a file
# Inputs:
#   -- SOURCE_DIR
# Outputs:
#   -- GEOS_VERSION_FULL
#   -- GEOS_VERSION_LIST
#   -- GEOS_VERSION_MAJOR
#   -- GEOS_VERSION_MINOR
#   -- GEOS_VERSION_PATCH
macro(geosx_get_file_version)
  file ( STRINGS "${SOURCE_DIR}/VERSION" GEOS_VERSION_FULL )
  string( REGEX REPLACE "VERSION_ID = v" "" GEOS_VERSION_FULL "${GEOS_VERSION_FULL}" )
  string( REPLACE "." ";" GEOS_VERSION_LIST ${GEOS_VERSION_FULL} )
  list( GET GEOS_VERSION_LIST  0 GEOS_VERSION_MAJOR )
  list( GET GEOS_VERSION_LIST  1 GEOS_VERSION_MINOR )
  list( GET GEOS_VERSION_LIST  2 GEOS_VERSION_PATCH )
endmacro()

# Get GEOSX development version from git
# Inputs:
#   -- SOURCE_DIR
#   -- GIT_FOUND
#   -- GIT_EXECUTABLE (only when GIT_FOUND=TRUE )
# Outputs (if GIT_FOUND=TRUE and inside git repo):
#   -- GEOS_GIT_BRANCH
#   -- GEOS_GIT_HASH
#   -- GEOS_GIT_TAG
macro(geosx_get_git_version)
  if( GIT_FOUND )
    # Use BLT Git macros for convenience
    include( ${SOURCE_DIR}/cmake/blt/cmake/BLTGitMacros.cmake )
    blt_is_git_repo( OUTPUT_STATE is_git_repo
                     SOURCE_DIR ${SOURCE_DIR} )
    if( is_git_repo )
      blt_git_branch( BRANCH_NAME GEOS_GIT_BRANCH
                      RETURN_CODE _git_rc
                      SOURCE_DIR ${SOURCE_DIR} )
      blt_git_hashcode( HASHCODE GEOS_GIT_HASH
                        RETURN_CODE _git_rc
                        SOURCE_DIR ${SOURCE_DIR} )
      blt_git_tag( OUTPUT_TAG GEOS_GIT_TAG
                   RETURN_CODE _git_rc
                   SOURCE_DIR ${SOURCE_DIR} )
    endif()
  endif()
endmacro()

# Both at config time and build time, we generate updated version info
geosx_get_file_version()
geosx_get_git_version()

# Set paths for version header
set( ver_component mainInterface )
set( ver_filename GeosxVersion.hpp )
set( ver_in_file ${SOURCE_DIR}/coreComponents/${ver_component}/${ver_filename}.in )
set( ver_tmp_file ${CMAKE_BINARY_DIR}/${ver_filename} )
set( ver_out_path ${CMAKE_BINARY_DIR}/include/${ver_component} )
set( ver_out_file ${ver_out_path}/${ver_filename} )

# Generate version header, but in a temporary location to avoid rebuilds every time
configure_file( ${ver_in_file} ${ver_tmp_file} )

if( NOT DEFINED CMAKE_SCRIPT_MODE_FILE )
  # At config time copy version header into place so it can be picked up by IDEs, etc.
  # Also make it installable in case another project wants to include it.
  file( COPY ${ver_tmp_file} DESTINATION ${ver_out_path} )
  install( FILES ${ver_out_file}
           DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${ver_component} )

  # Define a build-time target that will call this script again and:
  # - Regenerate the temporary version header
  # - Copy the header into place whenever it changes, to trigger partial rebuild
  add_custom_target( generate_version
                     ALL
                     DEPENDS ${ver_in_file}
                     BYPRODUCTS ${ver_out_file} ${ver_tmp_file}
                     COMMENT "Generating version information"
                     COMMAND ${CMAKE_COMMAND}
                     -D SOURCE_DIR=${SOURCE_DIR}
                     -D GIT_FOUND=${GIT_FOUND}
                     -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                     -P "${CMAKE_CURRENT_LIST_FILE}"
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different ${ver_tmp_file} ${ver_out_file} )
endif()
