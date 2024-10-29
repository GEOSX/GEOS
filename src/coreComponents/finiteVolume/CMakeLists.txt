#
# Specify all headers
#
set( finiteVolume_headers
     BoundaryStencil.hpp
     CellElementStencilMPFA.hpp
     CellElementStencilTPFA.hpp
     SurfaceElementStencil.hpp
     EmbeddedSurfaceToCellStencil.hpp
     FaceElementToCellStencil.hpp
     FiniteVolumeManager.hpp
     FluxApproximationBase.hpp
     ProjectionEDFMHelper.hpp
     HybridMimeticDiscretization.hpp
     MimeticInnerProductDispatch.hpp
     StencilBase.hpp
     TwoPointFluxApproximation.hpp
     mimeticInnerProducts/BdVLMInnerProduct.hpp
     mimeticInnerProducts/MimeticInnerProductBase.hpp
     mimeticInnerProducts/MimeticInnerProductHelpers.hpp
     mimeticInnerProducts/QuasiRTInnerProduct.hpp
     mimeticInnerProducts/QuasiTPFAInnerProduct.hpp
     mimeticInnerProducts/SimpleInnerProduct.hpp
     mimeticInnerProducts/TPFAInnerProduct.hpp
   )
#
# Specify all sources
#
set( finiteVolume_sources
     BoundaryStencil.cpp
     CellElementStencilMPFA.cpp
     CellElementStencilTPFA.cpp
     SurfaceElementStencil.cpp
     FaceElementToCellStencil.cpp
     EmbeddedSurfaceToCellStencil.cpp
     FiniteVolumeManager.cpp
     FluxApproximationBase.cpp
     TwoPointFluxApproximation.cpp
     ProjectionEDFMHelper.cpp
     HybridMimeticDiscretization.cpp
   )

set( dependencyList ${parallelDeps} mesh codingUtilities fieldSpecification )

blt_add_library( NAME       finiteVolume
                 SOURCES    ${finiteVolume_sources}
                 HEADERS    ${finiteVolume_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
               )

target_include_directories( finiteVolume PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

