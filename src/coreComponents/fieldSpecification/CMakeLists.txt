#
# Specify all headers
#
set( fieldSpecification_headers
     DirichletBoundaryCondition.hpp
     EquilibriumInitialCondition.hpp
     FieldSpecificationBase.hpp
     FieldSpecificationManager.hpp
     SourceFluxBoundaryCondition.hpp
     TractionBoundaryCondition.hpp
     AquiferBoundaryCondition.hpp
     PerfectlyMatchedLayer.hpp
   )

#
# Specify all sources
#
set( fieldSpecification_sources
     DirichletBoundaryCondition.cpp
     EquilibriumInitialCondition.cpp
     FieldSpecificationBase.cpp
     FieldSpecificationManager.cpp
     SourceFluxBoundaryCondition.cpp
     TractionBoundaryCondition.cpp
     AquiferBoundaryCondition.cpp
     PerfectlyMatchedLayer.cpp
   )

set( dependencyList ${parallelDeps} mesh )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       fieldSpecification
                 SOURCES    ${fieldSpecification_sources}
                 HEADERS    ${fieldSpecification_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( fieldSpecification PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS fieldSpecification LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
