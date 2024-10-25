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

#include "ElementRegionBase.hpp"

#include "common/TimingMacros.hpp"


namespace geos
{
using namespace dataRepository;


ElementRegionBase::ElementRegionBase( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_materialList(),
  m_meshBody()
{

  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  this->registerGroup( viewKeyStruct::elementSubRegions() );

  registerWrapper( viewKeyStruct::materialListString(), &m_materialList ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of materials present in this region" );

  registerWrapper( viewKeyStruct::meshBodyString(), &m_meshBody ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "" ).
    setDescription( "Mesh body that contains this region" );

}


ElementRegionBase::~ElementRegionBase()
{}


string ElementRegionBase::verifyMeshBodyName( Group const & meshBodies,
                                              string const & meshBodyBlockName )
{
  string meshBodyName = meshBodyBlockName;
  localIndex const numberOfMeshBodies = meshBodies.numSubGroups();

  GEOS_THROW_IF( numberOfMeshBodies == 0,
                 "No MeshBodies found in this problem, please check if correct input file is provided", InputError );

  if( numberOfMeshBodies == 1 )
  {
    string const & onlyMeshBodyName = meshBodies.getGroup( 0 ).getName();

    if( meshBodyName=="" )
    {
      meshBodyName = onlyMeshBodyName;
    }
    GEOS_THROW_IF_NE_MSG( onlyMeshBodyName,
                          meshBodyName,
                          "MeshBody specified does not match MeshBody in hierarchy.",
                          InputError );
  }
  else
  {
    bool meshBodyFound = false;
    meshBodies.forSubGroups( [&] ( Group const & meshBody )
    {
      if( meshBody.getName()==meshBodyName )
      {
        meshBodyFound = true;
      }
    } );
    GEOS_THROW_IF( !meshBodyFound,
                   "There are multiple MeshBodies in this problem, but the "
                   "specified MeshBody name "<<meshBodyName<<" was not found",
                   InputError );
  }

  return meshBodyName;
}


}
