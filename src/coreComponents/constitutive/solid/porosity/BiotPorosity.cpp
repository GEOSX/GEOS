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
 * @file BiotPorosity.cpp
 */

#include "BiotPorosity.hpp"
#include "PorosityFields.hpp"
#include "constitutive/solid/SolidBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


BiotPorosity::BiotPorosity( string const & name, Group * const parent ):
  PorosityBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultGrainBulkModulusString(), &m_defaultGrainBulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Grain bulk modulus" );

  registerWrapper( viewKeyStruct::defaultThermalExpansionCoefficientString(), &m_defaultThermalExpansionCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default thermal expansion coefficient" );

  registerWrapper( viewKeyStruct::useUniaxialFixedStressString(), &m_useUniaxialFixedStress ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag enabling uniaxial approximation in fixed stress update" );

  registerField( fields::porosity::biotCoefficient{}, &m_biotCoefficient ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Biot coefficient" );

  registerField( fields::porosity::grainBulkModulus{}, &m_grainBulkModulus ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Grain Bulk modulus." );

  registerField( fields::porosity::thermalExpansionCoefficient{}, &m_thermalExpansionCoefficient );

  registerField( fields::porosity::meanTotalStressIncrement_k{}, &m_meanTotalStressIncrement_k );

  registerField( fields::porosity::averageMeanTotalStressIncrement_k{}, &m_averageMeanTotalStressIncrement_k );

  registerWrapper( viewKeyStruct::solidBulkModulusString(), &m_bulkModulus ).
    setApplyDefaultValue( 1e-6 ).
    setDescription( "Solid bulk modulus" );

  registerWrapper( viewKeyStruct::solidShearModulusString(), &m_shearModulus ).
    setApplyDefaultValue( 1e-6 ).
    setDescription( "Solid shear modulus" );
}

void BiotPorosity::allocateConstitutiveData( dataRepository::Group & parent,
                                             localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_meanTotalStressIncrement_k.resize( 0, numConstitutivePointsPerParentIndex );
}

void BiotPorosity::postInputInitialization()
{
  PorosityBase::postInputInitialization();

  getWrapper< array1d< real64 > >( fields::porosity::thermalExpansionCoefficient::key() ).
    setApplyDefaultValue( m_defaultThermalExpansionCoefficient );

  // set results as array default values
  getWrapper< array1d< real64 > >( fields::porosity::grainBulkModulus::key() ).
    setApplyDefaultValue( m_defaultGrainBulkModulus );
}

void BiotPorosity::initializeState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView1d< real64 const > referencePorosity = m_referencePorosity;
  arrayView2d< real64 >             newPorosity = m_newPorosity;
  arrayView2d< real64 >             porosity_n  = m_porosity_n;
  arrayView2d< real64 >         initialPorosity = m_initialPorosity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      newPorosity[k][q]     = referencePorosity[k];
      porosity_n[k][q]      = referencePorosity[k];
      initialPorosity[k][q] = referencePorosity[k];
    }
  } );
}

void BiotPorosity::saveConvergedState() const
{
  PorosityBase::saveConvergedState();
  m_meanTotalStressIncrement_k.zero();
  m_averageMeanTotalStressIncrement_k.zero();
}

void BiotPorosity::ignoreConvergedState() const
{
  PorosityBase::ignoreConvergedState();
  m_meanTotalStressIncrement_k.zero();
  m_averageMeanTotalStressIncrement_k.zero();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, BiotPorosity, string const &, Group * const )
} /* namespace constitutive */
} /* namespace geos */
