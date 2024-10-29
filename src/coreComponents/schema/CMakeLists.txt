

#
# Specify all headers
#
set( schema_headers
     schemaUtilities.hpp )

#
# Specify all sources
#
set( schema_sources
     schemaUtilities.cpp )

set( dependencyList ${parallelDeps} dataRepository )

blt_add_library( NAME       schema
                 SOURCES    ${schema_sources}
                 HEADERS    ${schema_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
               )


target_include_directories( schema PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

