/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
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

namespace geosx
{

using namespace dataRepository;

FluxApproximationBase::FluxApproximationBase( string const & name, Group * const parent )
  : Group( name, parent ),
  m_lengthScale( 1.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::fieldNameString(), &m_fieldName ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Name of primary solution field" );

  registerWrapper( viewKeyStruct::coeffNameString(), &m_coeffName ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Name of coefficient field" );

  registerWrapper( viewKeyStruct::targetRegionsString(), &m_targetRegions ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "List of regions to build the stencil for" );

  registerWrapper( viewKeyStruct::areaRelativeToleranceString(), &m_areaRelTol ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0e-8 ).
    setDescription( "Relative tolerance for area calculations." );
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

        Group & stencilParentGroup = mesh.registerGroup( groupKeyStruct::stencilMeshGroupString() );
        Group & stencilGroup = stencilParentGroup.registerGroup( getName() );

        registerCellStencil( stencilGroup );

        registerFractureStencil( stencilGroup );
      }
    } );
  } );

  domain.forMeshBodies( [&]( MeshBody & meshBody )
  {
    meshBody.forMeshLevels( [&]( MeshLevel & mesh )
    {
      if( mesh.isShallowCopy() )
      {
        Group & parentMesh = mesh.getShallowParent();
        Group & parentStencilParentGroup = parentMesh.getGroup( groupKeyStruct::stencilMeshGroupString() );
        mesh.registerGroup( groupKeyStruct::stencilMeshGroupString(), &parentStencilParentGroup );
      }
    } );
  } );

  // domain.forMeshBodies( [&]( MeshBody & meshBody )
  // {
  //   meshBody.forMeshLevels( [&]( MeshLevel & mesh )
  //   {
  //     if( mesh.isShallowCopy() )
  //     {
  //       Group & parentMesh = mesh.getShallowParent();
  //       Group & parentStencilParentGroup = parentMesh.getGroup( groupKeyStruct::stencilMeshGroupString() );
  //       mesh.registerGroup( groupKeyStruct::stencilMeshGroupString(), &parentStencilParentGroup );
  //     }
  //   } );
  // } );

}

void FluxApproximationBase::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

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

        Group & stencilParentGroup = mesh.getGroup( groupKeyStruct::stencilMeshGroupString() );
        Group & stencilGroup = stencilParentGroup.getGroup( getName() );
        // For each face-based Dirichlet boundary condition on target field, create a boundary stencil
        // TODO: Apply() should take a MeshLevel directly
        fsManager.apply< FaceManager >( 0.0, // time = 0
                                        mesh,
                                        m_fieldName,
                                        [&] ( FieldSpecificationBase const &,
                                              string const & setName,
                                              SortedArrayView< localIndex const > const &,
                                              FaceManager const &,
                                              string const & )
        {
          registerBoundaryStencil( stencilGroup, setName );
        } );

        // For each aquifer boundary condition, create a boundary stencil
        fsManager.apply< FaceManager,
                         AquiferBoundaryCondition >( 0.0, // time = 0
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

        // Compute the main cell-based stencil
        computeCellStencil( mesh );

        // For each face-based boundary condition on target field, compute the boundary stencil weights
        fsManager.apply< FaceManager >( 0.0,
                                        mesh,
                                        m_fieldName,
                                        [&] ( FieldSpecificationBase const &,
                                              string const & setName,
                                              SortedArrayView< localIndex const > const & faceSet,
                                              FaceManager const &,
                                              string const & )
        {
          computeBoundaryStencil( mesh, setName, faceSet );
        } );

        // Compute the aquifer stencil weights
        computeAquiferStencil( domain, mesh );
      }
    } );
  } );
}

void FluxApproximationBase::setFieldName( string const & name )
{
  m_fieldName = name;
}

void FluxApproximationBase::setCoeffName( string const & name )
{
  m_coeffName = name;
}


} //namespace geosx
