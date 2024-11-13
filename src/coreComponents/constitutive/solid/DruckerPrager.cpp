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
 *  @file DruckerPrager.cpp
 */

#include "DruckerPrager.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

DruckerPrager::DruckerPrager( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_defaultFrictionAngle(),
  m_defaultDilationAngle(),
  m_defaultCohesion(),
  m_defaultHardening(),
  m_friction(),
  m_dilation(),
  m_hardening(),
  m_newCohesion(),
  m_oldCohesion()
{
  // register default values

  registerWrapper( viewKeyStruct::defaultFrictionAngleString(), &m_defaultFrictionAngle ).
    setApplyDefaultValue( 30.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Friction angle (degrees)" );

  registerWrapper( viewKeyStruct::defaultDilationAngleString(), &m_defaultDilationAngle ).
    setApplyDefaultValue( 30.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Dilation angle (degrees)" );

  registerWrapper( viewKeyStruct::defaultHardeningString(), &m_defaultHardening ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Cohesion hardening/softening rate" );

  registerWrapper( viewKeyStruct::defaultCohesionString(), &m_defaultCohesion ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Initial cohesion" );

  // register fields

  registerWrapper( viewKeyStruct::frictionString(), &m_friction ).
    setApplyDefaultValue( -1 ).
    setDescription( "Yield surface slope" );

  registerWrapper( viewKeyStruct::dilationString(), &m_dilation ).
    setApplyDefaultValue( -1 ).
    setDescription( "Plastic potential slope" );

  registerWrapper( viewKeyStruct::hardeningString(), &m_hardening ).
    setApplyDefaultValue( -1 ).
    setDescription( "Hardening rate" );

  registerWrapper( viewKeyStruct::newCohesionString(), &m_newCohesion ).
    setApplyDefaultValue( -1 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_3 ).
    setDescription( "New cohesion state" );

  registerWrapper( viewKeyStruct::oldCohesionString(), &m_oldCohesion ).
    setApplyDefaultValue( -1 ).
    setDescription( "Old cohesion state" );
}


DruckerPrager::~DruckerPrager()
{}


void DruckerPrager::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  m_newCohesion.resize( 0, numConstitutivePointsPerParentIndex );
  m_oldCohesion.resize( 0, numConstitutivePointsPerParentIndex );

  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


void DruckerPrager::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  GEOS_THROW_IF( m_defaultCohesion < 0,
                 getFullName() << ": Negative cohesion value detected", InputError );
  GEOS_THROW_IF( m_defaultFrictionAngle < 0,
                 getFullName() << ": Negative friction angle detected", InputError );
  GEOS_THROW_IF( m_defaultDilationAngle < 0,
                 getFullName() << ": Negative dilation angle detected", InputError );
  GEOS_THROW_IF( m_defaultFrictionAngle < m_defaultDilationAngle,
                 getFullName() << ": Dilation angle should not exceed friction angle", InputError );

  // convert from Mohr-Coulomb constants to Drucker-Prager constants, assuming DP
  // passes through the triaxial compression corners of the MC surface.
  // see Borja (2013) p. 75

  real64 phi = m_defaultFrictionAngle * M_PI / 180;
  real64 psi = m_defaultDilationAngle * M_PI / 180;

  real64 C = 6 * m_defaultCohesion * cos( phi ) / ( 3 - sin( phi ) );
  real64 F = 6 * sin( phi ) / ( 3 - sin( phi ) );
  real64 D = 6 * sin( psi ) / ( 3 - sin( psi ) );

  // set results as array default values

  getWrapper< array2d< real64 > >( viewKeyStruct::oldCohesionString() ).
    setApplyDefaultValue( C );

  getWrapper< array2d< real64 > >( viewKeyStruct::newCohesionString() ).
    setApplyDefaultValue( C );

  getWrapper< array1d< real64 > >( viewKeyStruct::dilationString() ).
    setApplyDefaultValue( D );

  getWrapper< array1d< real64 > >( viewKeyStruct::frictionString() ).
    setApplyDefaultValue( F );

  getWrapper< array1d< real64 > >( viewKeyStruct::hardeningString() ).
    setApplyDefaultValue( m_defaultHardening );
}


void DruckerPrager::saveConvergedState() const
{
  SolidBase::saveConvergedState(); // TODO: not ideal, as we have separate loops for base and derived data

  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView2d< real64 const > newCohesion = m_newCohesion;
  arrayView2d< real64 > oldCohesion = m_oldCohesion;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      oldCohesion( k, q ) = newCohesion( k, q );
    }
  } );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, DruckerPrager, std::string const &, Group * const )
}
} /* namespace geos */
