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
 * @file PlanarGeometricObject.cpp
 */

#include "PlanarGeometricObject.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
using namespace dataRepository;

PlanarGeometricObject::PlanarGeometricObject( const string & name, Group * const parent ):
  SimpleGeometricObjectBase( name, parent ),
  m_normal{ 0.0, 0.0, 1.0 },
  m_lengthVector{ 0.0, 0.0, 0.0 },
  m_widthVector{ 0.0, 0.0, 0.0 }
{
  registerWrapper( viewKeyStruct::normalString(), &m_normal ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)" );

  registerWrapper( viewKeyStruct::mLengthVectorString(), &m_lengthVector ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tangent vector defining the orthonormal basis along with the normal." );

  registerWrapper( viewKeyStruct::mWidthVectorString(), &m_widthVector ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tangent vector defining the orthonormal basis along with the normal." );
}

PlanarGeometricObject::~PlanarGeometricObject()
{}

//REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, PlanarGeometricObject, string const &, Group * const )

} /* namespace geos */
