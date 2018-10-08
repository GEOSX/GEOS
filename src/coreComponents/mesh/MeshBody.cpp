/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * MeshBody.cpp
 *
 *  Created on: Sep 13, 2017
 *      Author: settgast
 */

#include "MeshBody.hpp"
#include "MeshLevel.hpp"

namespace geosx
{
using namespace dataRepository;

MeshBody::MeshBody( string const & name,
                    ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  RegisterViewWrapper<integer>( viewKeys.meshLevels );
}

MeshBody::MeshBody( string const & name,
                    ManagedGroup * const parent,
                    string const & meshBodyGeneratorType) :
    MeshBody( name, parent)
{
}

MeshBody::~MeshBody()
{
  // TODO Auto-generated destructor stub
}



MeshLevel * MeshBody::CreateMeshLevel( integer const newLevel )
{
  return this->RegisterGroup<MeshLevel>( "Level0" );
}

MeshBody::CatalogInterface::CatalogType& MeshBody::GetCatalog()
{
  static MeshBody::CatalogInterface::CatalogType catalog;
  return catalog;
}

void MeshBody::CreateChild( string const & childKey, string const & childName )
{
    if( childKey == "Mesh" ) {
    }
    else if( childKey == "Properties") {
        //TODO : Implement this case ---> next PR !
    }
    else
    {
        GEOS_ERROR(childKey << "is not a valid key that can be nested in " << CatalogName());
    }
}

} /* namespace geosx */
