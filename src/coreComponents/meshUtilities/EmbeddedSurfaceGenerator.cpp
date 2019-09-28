/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file EmbeddedSurfaceGenerator.cpp
 */

#include "EmbeddedSurfaceGenerator.hpp"


namespace geosx
{

#include "mpiCommunications/CommunicationTools.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "EmbeddedSurfaceRegion.hpp"
#include "EmbeddedSurfaceSubRegion.hpp"

namespace geosx
{
using namespace dataRepository;

EmbeddedSurfaceGenerator::EmbeddedSurfaceGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent ),
  m_meshBodyName(""),
  m_numElems(0),
  m_numNodesPerElem(4),
  m_numNodes(0),
  m_nDims(3)
{
  registerWrapper(keys::nVector, &m_nVector, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("normal unit vector to the fracture plane.");

  registerWrapper(keys::planeCenter, &m_planeCenter, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("coordinates of the center of the plane defining the surface.");

  registerWrapper(keys::meshBodyName, &m_meshBodyName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("");
}

EmbeddedSurfaceGenerator::~EmbeddedSurfaceGenerator()
{
  // TODO Auto-generated destructor stub
}

void EmbeddedSurfaceGenerator::PostProcessInput()
{
  // GEOS_ERROR_IF( m_nVector.empty(),
     //              "Cannot be empy " << getName() );

  //GEOS_ERROR_IF( m_planeCenter.empty(),
    //               "Cannot be empty " << getName() );

  GEOS_ERROR_IF( m_meshBodyName.empty(),
                 "meshName cannot be empty " << getName() );
 }

Group * EmbeddedSurfaceGenerator::CreateChild( string const & GEOSX_UNUSED_ARG( childKey ), string const & GEOSX_UNUSED_ARG( childName ) )
{
  // does not have any children (at least for now)
  return nullptr;
}

void EmbeddedSurfaceGenerator::GenerateMesh( DomainPartition * const domain )
{
  /*
   * Here we need to:
   * 1. Find intersections between the fracture plane and the background mesh.
   * 2. Generate EmbeddedSurface elements and populate them
   */



}


REGISTER_CATALOG_ENTRY( MeshGeneratorBase,
                        EmbeddedSurfaceGenerator,
                        std::string const &, dataRepository::Group* const )

} /* namespace geosx */
