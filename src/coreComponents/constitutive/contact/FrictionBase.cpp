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
 * @file FrictionBase.cpp
 */

#include "FrictionBase.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

FrictionBase::FrictionBase( string const & name,
                            Group * const parent ):
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::displacementJumpThresholdString(), &m_displacementJumpThreshold ).
    setApplyDefaultValue( std::numeric_limits< real64 >::epsilon() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "A threshold valued to determine whether a fracture is open or not." );
}

FrictionBase::~FrictionBase()
{}

FrictionBaseUpdates FrictionBase::createKernelWrapper() const
{
  return FrictionBaseUpdates( m_displacementJumpThreshold );
}

} /* end namespace constitutive */

} /* end namespace geos */
