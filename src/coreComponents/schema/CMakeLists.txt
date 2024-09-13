

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

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       schema
                 SOURCES    ${schema_sources}
                 HEADERS    ${schema_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )


target_include_directories( schema PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

install( TARGETS schema LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
