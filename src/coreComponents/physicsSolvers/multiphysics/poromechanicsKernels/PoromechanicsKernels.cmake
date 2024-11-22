set( kernelPath "coreComponents/physicsSolvers/multiphysics/poromechanicsKernels" )

set( SinglePhasePoromechanicsPolicy "geos::parallelDevicePolicy< ${GEOS_BLOCK_SIZE} >" )
set( SinglePhasePoromechanicsEFEMPolicy "geos::parallelDevicePolicy< ${GEOS_BLOCK_SIZE} >" )
set( SinglePhasePoromechanicsDamagePolicy "geos::parallelDevicePolicy< ${GEOS_BLOCK_SIZE} >" )
set( MultiphasePoromechanicsPolicy "geos::parallelDevicePolicy< ${GEOS_BLOCK_SIZE} >" )
set( ThermalMultiphasePoromechanicsPolicy "geos::parallelDevicePolicy< ${GEOS_BLOCK_SIZE} >" )
set( ThermalSinglePhasePoromechanicsPolicy "geos::parallelDevicePolicy< ${GEOS_BLOCK_SIZE} >" )
set( ThermalSinglePhasePoromechanicsEFEMPolicy "geos::parallelDevicePolicy< ${GEOS_BLOCK_SIZE} >" )


configure_file( ${CMAKE_SOURCE_DIR}/${kernelPath}/policies.hpp.in
                ${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/policies.hpp )

set( kernelNames PoromechanicsKernels ThermoPoromechanicsKernels )
set( damageKernelNames PoromechanicsDamageKernels )
set( subregionList CellElementSubRegion )
set( porousSolidDispatch PorousSolid<DruckerPragerExtended>
                         PorousSolid<ModifiedCamClay>
                         PorousSolid<DelftEgg>
                         PorousSolid<DruckerPrager>
                         PorousSolid<ElasticIsotropic>
                         PorousSolid<ElasticTransverseIsotropic>
                         PorousSolid<ElasticIsotropicPressureDependent>
                         PorousSolid<ElasticOrthotropic>
                         PorousSolid<DamageSpectral<ElasticIsotropic>>
                         PorousSolid<DamageVolDev<ElasticIsotropic>>
                         PorousSolid<Damage<ElasticIsotropic>> 
                         PorousSolid<DuvautLionsSolid<DruckerPrager>>
                         PorousSolid<DuvautLionsSolid<DruckerPragerExtended>>
                         PorousSolid<DuvautLionsSolid<ModifiedCamClay>> )

set( porousDamageSolidDispatch PorousDamageSolid<Damage<ElasticIsotropic>>
                               PorousDamageSolid<DamageSpectral<ElasticIsotropic>>
                               PorousDamageSolid<DamageVolDev<ElasticIsotropic>> )

set( finiteElementDispatch H1_Hexahedron_Lagrange1_GaussLegendre2
                           H1_Wedge_Lagrange1_Gauss6
                           H1_Tetrahedron_Lagrange1_Gauss1
                           H1_Pyramid_Lagrange1_Gauss5
                           H1_Tetrahedron_VEM_Gauss1
                           H1_Prism5_VEM_Gauss1
                           H1_Prism6_VEM_Gauss1
                           H1_Prism7_VEM_Gauss1
                           H1_Prism8_VEM_Gauss1
                           H1_Prism9_VEM_Gauss1
                           H1_Prism10_VEM_Gauss1 )

if ( NOT ${ENABLE_HIP} )
  list(APPEND finiteElementDispatch
              H1_Hexahedron_VEM_Gauss1
              H1_Wedge_VEM_Gauss1
              H1_Prism11_VEM_Gauss1 )
endif( )


  foreach( KERNELNAME ${kernelNames} )
    foreach( SUBREGION_TYPE  ${subregionList} )
      foreach( CONSTITUTIVE_TYPE ${porousSolidDispatch} )
        foreach( FE_TYPE ${finiteElementDispatch} )

        set( filename "${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/${KERNELNAME}_${SUBREGION_TYPE}_${CONSTITUTIVE_TYPE}_${FE_TYPE}.cpp" )
        string(REPLACE "<" "-" filename ${filename})
        string(REPLACE ">" "-" filename ${filename})
        string(REPLACE "," "-" filename ${filename})
        string(REPLACE " " "" filename ${filename})
        message( " -- Generating file: ${filename}")
        configure_file( ${CMAKE_SOURCE_DIR}/${kernelPath}/${KERNELNAME}.cpp.template
                        ${filename} )

          list( APPEND physicsSolvers_sources ${filename} )
        endforeach()
      endforeach()
    endforeach()
  endforeach()

  foreach( KERNELNAME ${damageKernelNames} )
    foreach( SUBREGION_TYPE  ${subregionList} )
      foreach( CONSTITUTIVE_TYPE ${porousDamageSolidDispatch} )
        foreach( FE_TYPE ${finiteElementDispatch} )

        set( filename "${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/${KERNELNAME}_${SUBREGION_TYPE}_${CONSTITUTIVE_TYPE}_${FE_TYPE}.cpp" )
        string(REPLACE "<" "-" filename ${filename})
        string(REPLACE ">" "-" filename ${filename})
        string(REPLACE "," "-" filename ${filename})
        string(REPLACE " " "" filename ${filename})
        message( " -- Generating file: ${filename}")
        configure_file( ${CMAKE_SOURCE_DIR}/${kernelPath}/PoromechanicsDamageKernels.cpp.template
                        ${filename} )

          list( APPEND physicsSolvers_sources ${filename} )
        endforeach()
      endforeach()
    endforeach()
  endforeach()

set( kernelNames PoromechanicsEFEMKernels )
set( subregionList CellElementSubRegion )
set( porousSolidDispatch PorousSolid<ElasticIsotropic> )
set( finiteElementDispatch H1_Hexahedron_Lagrange1_GaussLegendre2
                           H1_Wedge_Lagrange1_Gauss6
                           H1_Tetrahedron_Lagrange1_Gauss1
                           H1_Pyramid_Lagrange1_Gauss5
                           H1_Tetrahedron_VEM_Gauss1
                           H1_Prism5_VEM_Gauss1
                           H1_Prism6_VEM_Gauss1
                           H1_Prism7_VEM_Gauss1
                           H1_Prism8_VEM_Gauss1
                           H1_Prism9_VEM_Gauss1
                           H1_Prism10_VEM_Gauss1 )

if ( NOT ${ENABLE_HIP} )
  list(APPEND finiteElementDispatch
              H1_Hexahedron_VEM_Gauss1
              H1_Wedge_VEM_Gauss1
              H1_Prism11_VEM_Gauss1 )
endif( )
  
  foreach( KERNELNAME ${kernelNames} )
    foreach( SUBREGION_TYPE  ${subregionList} )
      foreach( CONSTITUTIVE_TYPE ${porousSolidDispatch} )
        foreach( FE_TYPE ${finiteElementDispatch} )
  
        set( filename "${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/${KERNELNAME}_${SUBREGION_TYPE}_${CONSTITUTIVE_TYPE}_${FE_TYPE}.cpp" )
        string(REPLACE "<" "-" filename ${filename})
        string(REPLACE ">" "-" filename ${filename})
        string(REPLACE "," "-" filename ${filename})
        string(REPLACE " " "" filename ${filename})
        message( " -- Generating file: ${filename}")
        configure_file( ${CMAKE_SOURCE_DIR}/${kernelPath}/PoromechanicsEFEMKernels.cpp.template
                          ${filename} )
  
        list( APPEND physicsSolvers_sources ${filename} )
        endforeach()
      endforeach()
    endforeach()
  endforeach()  

