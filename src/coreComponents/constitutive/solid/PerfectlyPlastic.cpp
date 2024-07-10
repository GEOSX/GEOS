/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file PerfectlyPlastic.cpp
 */

#include "PerfectlyPlastic.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

PerfectlyPlastic::PerfectlyPlastic( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_defaultYieldStress(),
  m_yieldStress()
{
  // register default values
  registerWrapper( viewKeyStruct::defaultYieldStressString(), &m_defaultYieldStress ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default yield stress" );

  // register fields
  registerWrapper( viewKeyStruct::yieldStressString(), &m_yieldStress ).
    setApplyDefaultValue( -1 ).
    setDescription( "Array of element yield stresses" );
}


PerfectlyPlastic::~PerfectlyPlastic()
{}


void PerfectlyPlastic::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


void PerfectlyPlastic::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  GEOS_THROW_IF( m_defaultYieldStress < 0.0, "Negative yield stress detected", InputError );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::yieldStressString() ).setApplyDefaultValue( m_defaultYieldStress );
}


void PerfectlyPlastic::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, PerfectlyPlastic, std::string const &, Group * const )
}
} /* namespace geos */
