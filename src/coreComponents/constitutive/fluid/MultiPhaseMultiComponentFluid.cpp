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
 * @file MultiPhaseMultiComponentFluid.cpp
 */
#include "MultiPhaseMultiComponentFluid.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

using namespace PVTProps;

namespace
{
template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH > class
  TwoPhaseCatalogNames {};

template<> class
  TwoPhaseCatalogNames< PVTProps::BrineCO2Density,
                        PVTProps::BrineViscosity,
                        PVTProps::SpanWagnerCO2Density,
                        PVTProps::FenghourCO2Viscosity,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrineFluid"; }
};
} // end namespace

// provide a definition for catalogName()
template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
string MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::catalogName()
{
  return TwoPhaseCatalogNames< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::name();
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::
MultiPhaseMultiComponentFluid( string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::phasePVTParaFilesString(), &m_phasePVTParaFiles ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Names of the files defining the parameters of the viscosity and density models" );

  registerWrapper( viewKeyStruct::flashModelParaFileString(), &m_flashModelParaFile ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Name of the file defining the parameters of the flash model" );
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
std::unique_ptr< ConstitutiveBase >
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  MultiPhaseMultiComponentFluid * const newConstitutiveRelation = dynamic_cast< MultiPhaseMultiComponentFluid * >(clone.get());
  newConstitutiveRelation->m_p1Index = this->m_p1Index;
  newConstitutiveRelation->m_p2Index = this->m_p2Index;

  newConstitutiveRelation->createPVTModels();

  return clone;
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
void MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  GEOSX_THROW_IF_NE_MSG( numFluidPhases(), 2,
                         "CO2BrineFluid model named " << getName() << ": The number of phases in this model should be equal to 2",
                         InputError );
  GEOSX_THROW_IF_NE_MSG( numFluidComponents(), 2,
                         "CO2BrineFluid model named " << getName() << ": The number of components in this model should be equal to 2",
                         InputError );
  GEOSX_THROW_IF_NE_MSG( m_phasePVTParaFiles.size(), 2,
                         "CO2BrineFluid model named " << getName() << ": The number of phasePVTParaFiles is not the same as the number of phases!",
                         InputError );

  // NOTE: for now, the names of the phases are still hardcoded here
  // Later, we could read them from the XML file and we would then have a general class here

  string const expectedWaterPhaseNames[] = { "Water", "water", "Liquid", "liquid" };
  m_p1Index = PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );

  string const expectedGasPhaseNames[] = { "CO2", "co2", "gas", "Gas" };
  m_p2Index = PVTFunctionHelpers::findName( m_phaseNames, expectedGasPhaseNames, viewKeyStruct::phaseNamesString() );

  createPVTModels();
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
void MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::createPVTModels()
{
  // 1) Create the viscosity and density models
  for( string const & filename : m_phasePVTParaFiles )
  {
    std::ifstream is( filename );
    string str;
    while( std::getline( is, str ) )
    {
      string_array const strs = stringutilities::tokenize( str, " " );

      if( strs[0] == "DensityFun" )
      {
        if( strs[1] == P1DENS::catalogName() )
        {
          m_p1Density = std::make_unique< P1DENS >( strs, m_componentNames, m_componentMolarWeight );
        }
        else if( strs[1] == P2DENS::catalogName() )
        {
          m_p2Density = std::make_unique< P2DENS >( strs, m_componentNames, m_componentMolarWeight );
        }
      }
      else if( strs[0] == "ViscosityFun" )
      {
        if( strs[1] == P1VISC::catalogName() )
        {
          m_p1Viscosity = std::make_unique< P1VISC >( strs, m_componentNames, m_componentMolarWeight );
        }
        else if( strs[1] == P2VISC::catalogName() )
        {
          m_p2Viscosity = std::make_unique< P2VISC >( strs, m_componentNames, m_componentMolarWeight );
        }
      }
      else
      {
        GEOSX_THROW( "CO2BrineFluid model named " << getName() << ": Invalid PVT function: " << strs[0] << ".", InputError );
      }
    }
    is.close();
  }

  GEOSX_THROW_IF( m_p1Density == nullptr, "CO2BrineFluid model named " << getName() << ": " << P1DENS::catalogName() << " model not found", InputError );
  GEOSX_THROW_IF( m_p2Density == nullptr, "CO2BrineFluid model named " << getName() << ": " << P2DENS::catalogName() << " model not found", InputError );
  GEOSX_THROW_IF( m_p1Viscosity == nullptr, "CO2BrineFluid model named " << getName() << ": " << P1VISC::catalogName() << " model not found", InputError );
  GEOSX_THROW_IF( m_p2Viscosity == nullptr, "CO2BrineFluid model named " << getName() << ": " << P2VISC::catalogName() << " model not found", InputError );

  // 2) Create the flash model
  {
    std::ifstream is( m_flashModelParaFile );
    string str;
    while( std::getline( is, str ) )
    {
      string_array const strs = stringutilities::tokenize( str, " " );
      if( strs[0] == "FlashModel" && strs[1] == FLASH::catalogName() )
      {
        m_flash = std::make_unique< FLASH >( strs, m_phaseNames, m_componentNames, m_componentMolarWeight );
      }
      else
      {
        GEOSX_THROW( "CO2BrineFluid model named " << getName() << ": Invalid flash model: " << strs[0] << ", " << strs[1] << ".", InputError );
      }
    }
    is.close();
  }

  GEOSX_THROW_IF( m_flash == nullptr,
                  "CO2BrineFluid model named " << getName() << ": " << FLASH::catalogName() << " not found",
                  InputError );
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
typename MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::KernelWrapper
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::createKernelWrapper() const
{
  return KernelWrapper( m_p1Index,
                        m_p2Index,
                        *m_p1Density,
                        *m_p1Viscosity,
                        *m_p2Density,
                        *m_p2Viscosity,
                        *m_flash,
                        m_componentMolarWeight.toViewConst(),
                        m_useMass,
                        { m_phaseFraction,
                          m_dPhaseFraction_dPressure,
                          m_dPhaseFraction_dTemperature,
                          m_dPhaseFraction_dGlobalCompFraction },
                        { m_phaseDensity,
                          m_dPhaseDensity_dPressure,
                          m_dPhaseDensity_dTemperature,
                          m_dPhaseDensity_dGlobalCompFraction },
                        { m_phaseMassDensity,
                          m_dPhaseMassDensity_dPressure,
                          m_dPhaseMassDensity_dTemperature,
                          m_dPhaseMassDensity_dGlobalCompFraction },
                        { m_phaseViscosity,
                          m_dPhaseViscosity_dPressure,
                          m_dPhaseViscosity_dTemperature,
                          m_dPhaseViscosity_dGlobalCompFraction },
                        { m_phaseCompFraction,
                          m_dPhaseCompFraction_dPressure,
                          m_dPhaseCompFraction_dTemperature,
                          m_dPhaseCompFraction_dGlobalCompFraction },
                        { m_totalDensity,
                          m_dTotalDensity_dPressure,
                          m_dTotalDensity_dTemperature,
                          m_dTotalDensity_dGlobalCompFraction } );
}

template< typename P1DENS, typename P1VISC, typename P2DENS, typename P2VISC, typename FLASH >
MultiPhaseMultiComponentFluid< P1DENS, P1VISC, P2DENS, P2VISC, FLASH >::KernelWrapper::
  KernelWrapper( localIndex const p1Index,
                 localIndex const p2Index,
                 P1DENS const & p1DensityWrapper,
                 P1VISC const & p1ViscosityWrapper,
                 P2DENS const & p2DensityWrapper,
                 P2VISC const & p2ViscosityWrapper,
                 FLASH const & flashWrapper,
                 arrayView1d< geosx::real64 const > const & componentMolarWeight,
                 bool useMass,
                 MultiFluidBase::KernelWrapper::PhasePropViews const & phaseFraction,
                 MultiFluidBase::KernelWrapper::PhasePropViews const & phaseDensity,
                 MultiFluidBase::KernelWrapper::PhasePropViews const & phaseMassDensity,
                 MultiFluidBase::KernelWrapper::PhasePropViews const & phaseViscosity,
                 MultiFluidBase::KernelWrapper::PhaseCompViews const & phaseCompFraction,
                 MultiFluidBase::KernelWrapper::FluidPropViews const & totalDensity )
  : MultiFluidBase::KernelWrapper( componentMolarWeight,
                                   useMass,
                                   phaseFraction,
                                   phaseDensity,
                                   phaseMassDensity,
                                   phaseViscosity,
                                   phaseCompFraction,
                                   totalDensity ),
  m_p1Index( p1Index ),
  m_p2Index( p2Index ),
  m_p1Density( p1DensityWrapper.createKernelWrapper() ),
  m_p1Viscosity( p1ViscosityWrapper.createKernelWrapper() ),
  m_p2Density( p2DensityWrapper.createKernelWrapper() ),
  m_p2Viscosity( p2ViscosityWrapper.createKernelWrapper() ),
  m_flash( flashWrapper.createKernelWrapper() )
{}

// explicit instantiation of the model template; unfortunately we can't use CO2BrineFluid alias for this
template class MultiPhaseMultiComponentFluid< PVTProps::BrineCO2Density,
                                              PVTProps::BrineViscosity,
                                              PVTProps::SpanWagnerCO2Density,
                                              PVTProps::FenghourCO2Viscosity,
                                              PVTProps::CO2Solubility >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrineFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geosx
