set(catalyst_headers
  CatalystActions.hpp
  ConduitCapsule.hpp
  GenericConduitCapsule.hpp
)

set(catalyst_sources
  CatalystActions.cpp
  ConduitCapsule.cpp
)

set(catalyst_template_sources
  GenericConduitCapsule.tpp
)

set(catalyst_dependencies
  catalyst::catalyst
)

list(APPEND ${catalyst_headers} ${catalyst_template_sources})

blt_add_library( NAME       fileIOCatalyst
                 SOURCES    ${catalyst_sources}
                 HEADERS    ${catalyst_headers}
                 DEPENDS_ON ${catalyst_dependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
               )

geosx_add_code_checks(PREFIX fileIOCatalyst )
