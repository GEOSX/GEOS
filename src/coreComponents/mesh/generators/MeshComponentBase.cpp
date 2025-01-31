/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "MeshComponentBase.hpp"

namespace geos
{

using namespace dataRepository;

MeshComponentBase::MeshComponentBase( const string & name,
                                      Group * const parent )
  : Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

MeshComponentBase::~MeshComponentBase()
{}

MeshComponentBase::CatalogInterface::CatalogType & MeshComponentBase::getCatalog()
{
  static CatalogInterface::CatalogType catalog;
  return catalog;
}

}
