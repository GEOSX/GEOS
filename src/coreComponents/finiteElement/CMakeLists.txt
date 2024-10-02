#
# Specify all headers
#
set( finiteElement_headers
     BilinearFormUtilities.hpp
     FiniteElementDiscretization.hpp
     FiniteElementDiscretizationManager.hpp
     FiniteElementDispatch.hpp
     Kinematics.h
     LinearFormUtilities.hpp
     PDEUtilities.hpp
     elementFormulations/FiniteElementBase.hpp
     elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp
     elementFormulations/H1_Pyramid_Lagrange1_Gauss5.hpp
     elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp
     elementFormulations/H1_Tetrahedron_Lagrange1_Gauss1.hpp
     elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp
     elementFormulations/H1_Wedge_Lagrange1_Gauss6.hpp
     elementFormulations/ConformingVirtualElementOrder1.hpp
     elementFormulations/ConformingVirtualElementOrder1_impl.hpp
     elementFormulations/LagrangeBasis1.hpp
     elementFormulations/LagrangeBasis2.hpp
     elementFormulations/LagrangeBasis3GL.hpp
     elementFormulations/LagrangeBasis5GL.hpp
     kernelInterface/InterfaceKernelBase.hpp
     kernelInterface/ImplicitKernelBase.hpp
     kernelInterface/KernelBase.hpp
     kernelInterface/SparsityKernelBase.hpp
   )

if( NOT ${ENABLE_HIP} )
  list(APPEND finiteElement_headers elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp )
endif()

#
# Specify all sources
#
set( finiteElement_sources
     FiniteElementDiscretization.cpp
     FiniteElementDiscretizationManager.cpp
   )

set( dependencyList ${parallelDeps} dataRepository )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       finiteElement
                 SOURCES    ${finiteElement_sources}
                 HEADERS    ${finiteElement_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( finiteElement PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS finiteElement LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )

