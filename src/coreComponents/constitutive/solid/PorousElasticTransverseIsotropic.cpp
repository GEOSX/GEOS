/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file PorousElasticTransverseIsotropic.cpp
 */

#include "PorousElasticTransverseIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{
PorousElasticTransverseIsotropic::PorousElasticTransverseIsotropic( string const & name,
                                                                    Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >( name, parent )
{}

PorousElasticTransverseIsotropic::~PorousElasticTransverseIsotropic() = default;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticTransverseIsotropic, string const &, Group * const )

}
} /* namespace geosx */
