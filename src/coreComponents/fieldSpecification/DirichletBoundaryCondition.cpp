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
 * @file DirichletBoundaryCondition.cpp
 */

#include "DirichletBoundaryCondition.hpp"

namespace geos
{
using namespace dataRepository;

DirichletBoundaryCondition::DirichletBoundaryCondition( string const & name, Group * const parent ):
  FieldSpecificationBase( name, parent )
{
  // TODO Auto-generated constructor stub

}

DirichletBoundaryCondition::~DirichletBoundaryCondition()
{
  // TODO Auto-generated destructor stub
}



REGISTER_CATALOG_ENTRY( FieldSpecificationBase, DirichletBoundaryCondition, string const &, Group * const )

} /* namespace geos */
