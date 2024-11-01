/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FluxApproximationBase.cpp
 *
 */

#include "FluxApproximationBase.hpp"

#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{

using namespace dataRepository;

FluxApproximationBase::FluxApproximationBase( string const & name, Group * const parent )
  : Group( name, parent ),
  m_lengthScale( 1.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::fieldNameString(), &m_fieldNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Name of primary solution field" ).setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::coeffNameString(), &m_coeffName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Name of coefficient field" );

  registerWrapper( viewKeyStruct::targetRegionsString(), &m_targetRegions ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "List of regions to build the stencil for" );

  registerWrapper( viewKeyStruct::areaRelativeToleranceString(), &m_areaRelTol ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0e-8 ).
    setDescription( "Relative tolerance for area calculations." );

  registerWrapper( viewKeyStruct::upwindingSchemeString(), &m_upwindingParams.upwindingScheme ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( UpwindingScheme::PPU ).
    setDescription( "Type of upwinding scheme. "
                    "Valid options:\n* " + EnumStrings< UpwindingScheme >::concat( "\n* " ) );

}

FluxApproximationBase::CatalogInterface::CatalogType &
FluxApproximationBase::getCatalog()
{
  static FluxApproximationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void FluxApproximationBase::initializePreSubGroups()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  domain.forMeshBodies( [&]( MeshBody & meshBody )
  {
    meshBody.forMeshLevels( [&]( MeshLevel & mesh )
    {
      // Proceed with regular procedure only if the MeshLevel is not a shallow copy
      if( !(mesh.isShallowCopy() ) )
      {
        // Group structure: mesh1/finiteVolumeStencils/myTPFA
        Group * stencilParentGroup = mesh.getGroupPointer( groupKeyStruct::stencilMeshGroupString() );
        // There can be more than one FluxApproximation object so we check if the the group has
        // already been registered.
        if( stencilParentGroup == nullptr )
        {
          stencilParentGroup = &(mesh.registerGroup( groupKeyStruct::stencilMeshGroupString() ));
        }

        Group & stencilGroup = stencilParentGroup->registerGroup( getName() );

        registerCellStencil( stencilGroup );

        registerFractureStencil( stencilGroup );
      }
      else
      {
        // There can be more than one FluxApproximation object so we check if the the group has
        // already been registered.
        if( !mesh.hasGroup( groupKeyStruct::stencilMeshGroupString() ) )
        {
          Group & parentMesh = mesh.getShallowParent();
          Group & parentStencilParentGroup = parentMesh.getGroup( groupKeyStruct::stencilMeshGroupString() );
          mesh.registerGroup( groupKeyStruct::stencilMeshGroupString(), &parentStencilParentGroup );
        }
      }
    } );
  } );
}

void FluxApproximationBase::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  domain.forMeshBodies( [&]( MeshBody & meshBody )
  {
    m_lengthScale = meshBody.getGlobalLengthScale();
    meshBody.forMeshLevels( [&]( MeshLevel & mesh )
    {
      if( !(mesh.isShallowCopy() ) )
      {
        // Group structure: mesh1/finiteVolumeStencils/myTPFA

        // Compute the main cell-based stencil
        computeCellStencil( mesh );

        // Compute the fracture related stencils (within the fracture itself,
        // but between the fracture and the matrix as well).
        computeFractureStencil( mesh );

        Group & stencilParentGroup = mesh.getGroup( groupKeyStruct::stencilMeshGroupString() );
        Group & stencilGroup = stencilParentGroup.getGroup( getName() );
        // For each face-based Dirichlet boundary condition on target field, create a boundary stencil
        // TODO: Apply() should take a MeshLevel directly
        for( auto const & fieldName : m_fieldNames )
        {

          fsManager.apply< FaceManager >( 0.0, // time = 0
                                          mesh,
                                          fieldName,
                                          [&] ( FieldSpecificationBase const &,
                                                string const & setName,
                                                SortedArrayView< localIndex const > const & faceSet,
                                                FaceManager const &,
                                                string const & )
          {
            if( !stencilGroup.hasWrapper( setName ) )
            {
              registerBoundaryStencil( stencilGroup, setName );
              computeBoundaryStencil( mesh, setName, faceSet );
            }
          } );
        }
        // For each aquifer boundary condition, create a boundary stencil
        fsManager.apply< FaceManager,
                         AquiferBoundaryCondition >( 0.0,   // time = 0
                                                     mesh,
                                                     AquiferBoundaryCondition::catalogName(),
                                                     [&] ( AquiferBoundaryCondition const &,
                                                           string const & setName,
                                                           SortedArrayView< localIndex const > const &,
                                                           FaceManager const &,
                                                           string const & )
        {
          registerAquiferStencil( stencilGroup, setName );
        } );
        // Compute the aquifer stencil weights
        computeAquiferStencil( domain, mesh );
      }
    } );
  } );
}

void FluxApproximationBase::addFieldName( string const & name )
{
  m_fieldNames.emplace_back( name );
}

void FluxApproximationBase::setCoeffName( string const & name )
{
  m_coeffName = name;
}


} //namespace geos
