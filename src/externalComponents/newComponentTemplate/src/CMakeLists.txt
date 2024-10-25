# Specify all headers
set( componentHeaders NewComponent.hpp )

# Specify all sources
set( componentSources NewComponent.cpp )

# blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)

unset( dependencyList )

# Specify dependencies
list( APPEND dependencyList geosx_core)

message("adding newComponentTemplate library")

blt_add_library( NAME             newComponentTemplate
                 SOURCES          ${componentSources}
                 HEADERS          ${componentHeaders}
                 DEPENDS_ON       ${dependencyList}
                 OBJECT           ${GEOS_BUILD_OBJ_LIBS} )

target_include_directories( newComponentTemplate PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
