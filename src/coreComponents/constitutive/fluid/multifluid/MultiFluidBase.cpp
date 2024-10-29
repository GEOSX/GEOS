/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiFluidBase.cpp
 */


#include "MultiFluidBase.hpp"
#include "MultiFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

MultiFluidBase::MultiFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent ),
  m_useMass( false ),
  m_checkPVTTablesRanges( 1 )
{
  // We make base inputs optional here, since derived classes may want to predefine/hardcode
  // components/phases. Models that do need these inputs should change input flags accordingly.

  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of component names" );

  registerWrapper( viewKeyStruct::componentMolarWeightString(), &m_componentMolarWeight ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component molar weights" );

  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::useMassString(), &m_useMass ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerField( fields::multifluid::phaseFraction{}, &m_phaseFraction.value );
  registerField( fields::multifluid::dPhaseFraction{}, &m_phaseFraction.derivs );

  registerField( fields::multifluid::phaseDensity{}, &m_phaseDensity.value );
  registerField( fields::multifluid::phaseDensity_n{}, &m_phaseDensity_n );
  registerField( fields::multifluid::dPhaseDensity{}, &m_phaseDensity.derivs );

  registerField( fields::multifluid::phaseMassDensity{}, &m_phaseMassDensity.value );
  registerField( fields::multifluid::dPhaseMassDensity{}, &m_phaseMassDensity.derivs );

  registerField( fields::multifluid::phaseViscosity{}, &m_phaseViscosity.value );
  registerField( fields::multifluid::dPhaseViscosity{}, &m_phaseViscosity.derivs );

  registerField( fields::multifluid::phaseEnthalpy{}, &m_phaseEnthalpy.value );
  registerField( fields::multifluid::phaseEnthalpy_n{}, &m_phaseEnthalpy_n );
  registerField( fields::multifluid::dPhaseEnthalpy{}, &m_phaseEnthalpy.derivs );

  registerField( fields::multifluid::phaseInternalEnergy{}, &m_phaseInternalEnergy.value );
  registerField( fields::multifluid::phaseInternalEnergy_n{}, &m_phaseInternalEnergy_n );
  registerField( fields::multifluid::dPhaseInternalEnergy{}, &m_phaseInternalEnergy.derivs );

  registerField( fields::multifluid::phaseCompFraction{}, &m_phaseCompFraction.value );
  registerField( fields::multifluid::phaseCompFraction_n{}, &m_phaseCompFraction_n );
  registerField( fields::multifluid::dPhaseCompFraction{}, &m_phaseCompFraction.derivs );

  registerField( fields::multifluid::totalDensity{}, &m_totalDensity.value );
  registerField( fields::multifluid::totalDensity_n{}, &m_totalDensity_n );
  registerField( fields::multifluid::dTotalDensity{}, &m_totalDensity.derivs );


  registerWrapper( viewKeyStruct::checkPVTTablesRangesString(), &m_checkPVTTablesRanges ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Enable (1) or disable (0) an error when the input pressure or temperature of the PVT tables is out of range." ).
    setDefaultValue( 1 );
}

void MultiFluidBase::resizeFields( localIndex const size, localIndex const numPts )
{
  integer const numPhase = numFluidPhases();
  integer const numComp = numFluidComponents();
  integer const numDof = numComp + 2;

  m_phaseFraction.value.resize( size, numPts, numPhase );
  m_phaseFraction.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseDensity.value.resize( size, numPts, numPhase );
  m_phaseDensity_n.resize( size, numPts, numPhase );
  m_phaseDensity.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseMassDensity.value.resize( size, numPts, numPhase );
  m_phaseMassDensity.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseViscosity.value.resize( size, numPts, numPhase );
  m_phaseViscosity.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseEnthalpy.value.resize( size, numPts, numPhase );
  m_phaseEnthalpy_n.resize( size, numPts, numPhase );
  m_phaseEnthalpy.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseInternalEnergy.value.resize( size, numPts, numPhase );
  m_phaseInternalEnergy_n.resize( size, numPts, numPhase );
  m_phaseInternalEnergy.derivs.resize( size, numPts, numPhase, numDof );

  m_phaseCompFraction.value.resize( size, numPts, numPhase, numComp );
  m_phaseCompFraction_n.resize( size, numPts, numPhase, numComp );
  m_phaseCompFraction.derivs.resize( size, numPts, numPhase, numComp, numDof );

  m_totalDensity.value.resize( size, numPts );
  m_totalDensity_n.resize( size, numPts );
  m_totalDensity.derivs.resize( size, numPts, numDof );
}

void MultiFluidBase::setLabels()
{
  getField< fields::multifluid::phaseFraction >().
    setDimLabels( 2, m_phaseNames );

  getField< fields::multifluid::phaseDensity >().
    setDimLabels( 2, m_phaseNames );

  getField< fields::multifluid::phaseMassDensity >().
    setDimLabels( 2, m_phaseNames );

  getField< fields::multifluid::phaseViscosity >().
    setDimLabels( 2, m_phaseNames );

  getField< fields::multifluid::phaseEnthalpy >().
    setDimLabels( 2, m_phaseNames );

  getField< fields::multifluid::phaseInternalEnergy >().
    setDimLabels( 2, m_phaseNames );

  getField< fields::multifluid::phaseCompFraction >().
    setDimLabels( 2, m_phaseNames ).
    setDimLabels( 3, m_componentNames );
}

void MultiFluidBase::allocateConstitutiveData( dataRepository::Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );
}

void MultiFluidBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  integer const numComp = numFluidComponents();
  integer const numPhase = numFluidPhases();

  GEOS_THROW_IF_LT_MSG( numComp, 1,
                        GEOS_FMT( "{}: invalid number of components", getFullName() ),
                        InputError );
  GEOS_THROW_IF_GT_MSG( numComp, MAX_NUM_COMPONENTS,
                        GEOS_FMT( "{}: invalid number of components", getFullName() ),
                        InputError );
  GEOS_THROW_IF_LT_MSG( numPhase, 1,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );
  GEOS_THROW_IF_GT_MSG( numPhase, MAX_NUM_PHASES,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );
  GEOS_THROW_IF_NE_MSG( m_componentMolarWeight.size(), numComp,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::componentMolarWeightString() ),
                        InputError );

  // call to correctly set member array tertiary sizes on the 'main' material object
  resizeFields( 0, 0 );

  // set labels on array wrappers for plottable fields
  setLabels();
}

void MultiFluidBase::initializeState() const
{
  // initialize the "old" variables
  saveConvergedState();
}

void MultiFluidBase::saveConvergedState() const
{
  localIndex const numElem = m_phaseMassDensity.value.size( 0 );
  localIndex const numGauss = m_phaseMassDensity.value.size( 1 );
  integer const numPhase = m_phaseMassDensity.value.size( 2 );
  integer const numComp = m_phaseCompFraction.value.size( 3 );

  FluidProp::ViewTypeConst const totalDensity = m_totalDensity.toViewConst();
  PhaseProp::ViewTypeConst const phaseDensity = m_phaseDensity.toViewConst();
  PhaseProp::ViewTypeConst const phaseEnthalpy = m_phaseEnthalpy.toViewConst();
  PhaseProp::ViewTypeConst const phaseInternalEnergy = m_phaseInternalEnergy.toViewConst();
  PhaseComp::ViewTypeConst const phaseCompFraction = m_phaseCompFraction.toViewConst();

  arrayView2d< real64, multifluid::USD_FLUID > const totalDensity_n = m_totalDensity_n.toView();
  arrayView3d< real64, multifluid::USD_PHASE > const phaseDensity_n = m_phaseDensity_n.toView();
  arrayView3d< real64, multifluid::USD_PHASE > const phaseEnthalpy_n = m_phaseEnthalpy_n.toView();
  arrayView3d< real64, multifluid::USD_PHASE > const phaseInternalEnergy_n = m_phaseInternalEnergy_n.toView();
  arrayView4d< real64, multifluid::USD_PHASE_COMP > const phaseCompFraction_n = m_phaseCompFraction_n.toView();

  forAll< parallelDevicePolicy<> >( numElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numGauss; ++q )
    {
      totalDensity_n[k][q] = totalDensity.value[k][q];
      for( integer ip = 0; ip < numPhase; ++ip )
      {
        phaseDensity_n[k][q][ip] = phaseDensity.value[k][q][ip];
        phaseEnthalpy_n[k][q][ip] = phaseEnthalpy.value[k][q][ip];
        phaseInternalEnergy_n[k][q][ip] = phaseInternalEnergy.value[k][q][ip];
        for( integer ic = 0; ic < numComp; ++ic )
        {
          phaseCompFraction_n[k][q][ip][ic] = phaseCompFraction.value[k][q][ip][ic];
        }
      }
    }
  } );


}


} // namespace constitutive

} // namespace geos
