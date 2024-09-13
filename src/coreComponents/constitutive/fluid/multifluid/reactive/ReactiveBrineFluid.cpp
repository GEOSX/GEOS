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
 * @file ReactiveBrineFluid.cpp
 */
#include "ReactiveBrineFluid.hpp"

#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "common/Units.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

using namespace PVTProps;

namespace
{
template< typename PHASE > class
  ReactiveBrineCatalogNames {};

template<> class
  ReactiveBrineCatalogNames< PhaseModel< PVTProps::WaterDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction > >
{
public:
  static string name() { return "ReactiveBrine"; }
};
template<> class
  ReactiveBrineCatalogNames< PhaseModel< PVTProps::WaterDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy > >
{
public:
  static string name() { return "ReactiveBrineThermal"; }
};

} // end namespace

// provide a definition for catalogName()
template< typename PHASE >
string ReactiveBrineFluid< PHASE > ::catalogName()
{
  return ReactiveBrineCatalogNames< PHASE > ::name();
}

template< typename PHASE >
ReactiveBrineFluid< PHASE > ::
ReactiveBrineFluid( string const & name, Group * const parent ):
  ReactiveMultiFluid( name, parent )
{
  registerWrapper( viewKeyStruct::phasePVTParaFilesString(), &m_phasePVTParaFiles ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Names of the files defining the parameters of the viscosity and density models" );

  // if this is a thermal model, we need to make sure that the arrays will be properly displayed and saved to restart
  if( isThermal() )
  {
    getField< fields::multifluid::phaseEnthalpy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );

    getField< fields::multifluid::phaseInternalEnergy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );
  }
}

template< typename PHASE >
bool ReactiveBrineFluid< PHASE > ::isThermal() const
{
  return ( PHASE::Enthalpy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() );
}


template< typename PHASE >
std::unique_ptr< ConstitutiveBase >
ReactiveBrineFluid< PHASE > ::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ReactiveMultiFluid::deliverClone( name, parent );

  ReactiveBrineFluid & newConstitutiveRelation = dynamicCast< ReactiveBrineFluid & >( *clone );

  newConstitutiveRelation.createPVTModels();

  return clone;
}

template< typename PHASE >
integer ReactiveBrineFluid< PHASE > ::getWaterPhaseIndex() const
{
  // There is only 1 phase
  return 0;
}


template< typename PHASE >
void ReactiveBrineFluid< PHASE > ::postInputInitialization()
{
  ReactiveMultiFluid::postInputInitialization();

  GEOS_THROW_IF_NE_MSG( numFluidPhases(), 1,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );
  GEOS_THROW_IF_NE_MSG( m_phasePVTParaFiles.size(), 1,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName() ),
                        InputError );

  createPVTModels();
}

template< typename PHASE >
void ReactiveBrineFluid< PHASE > ::createPVTModels()
{

  // TODO: get rid of these external files and move into XML, this is too error prone
  // For now, to support the legacy input, we read all the input parameters at once in the arrays below, and then we create the models
  array1d< array1d< string > > phase1InputParams;
  phase1InputParams.resize( 3 );

  // 1) Create the viscosity, density, enthalpy models
  for( string const & filename : m_phasePVTParaFiles )
  {
    std::ifstream is( filename );
    string str;
    while( std::getline( is, str ) )
    {
      array1d< string > const strs = stringutilities::tokenizeBySpaces< array1d >( str );

      if( !strs.empty() )
      {
        GEOS_THROW_IF( strs.size() < 2,
                       GEOS_FMT( "{}: missing PVT model in line '{}'", getFullName(), str ),
                       InputError );

        if( strs[0] == "DensityFun" )
        {
          if( strs[1] == PHASE::Density::catalogName() )
          {
            phase1InputParams[PHASE::InputParamOrder::DENSITY] = strs;
          }
        }
        else if( strs[0] == "ViscosityFun" )
        {
          if( strs[1] == PHASE::Viscosity::catalogName() )
          {
            phase1InputParams[PHASE::InputParamOrder::VISCOSITY] = strs;
          }
        }
        else if( strs[0] == "EnthalpyFun" )
        {
          if( strs[1] == PHASE::Enthalpy::catalogName() )
          {
            phase1InputParams[PHASE::InputParamOrder::ENTHALPY] = strs;
          }
        }
        else
        {
          GEOS_THROW( GEOS_FMT( "{}: invalid PVT function type '{}'", getFullName(), strs[0] ), InputError );
        }
      }
    }
    is.close();
  }

  // at this point, we have read the file and we check the consistency of non-thermal models
  GEOS_THROW_IF( phase1InputParams[PHASE::InputParamOrder::DENSITY].empty(),
                 GEOS_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE::Density::catalogName() ),
                 InputError );
  GEOS_THROW_IF( phase1InputParams[PHASE::InputParamOrder::VISCOSITY].empty(),
                 GEOS_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE::Viscosity::catalogName() ),
                 InputError );
  // we also detect any inconsistency arising in the enthalpy models
  GEOS_THROW_IF( phase1InputParams[PHASE::InputParamOrder::ENTHALPY].empty() &&
                 ( PHASE::Enthalpy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() ),
                 GEOS_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE::Enthalpy::catalogName() ),
                 InputError );

  // then, we are ready to instantiate the phase models
  m_phase = std::make_unique< PHASE >( getName() + "_phaseModel1", phase1InputParams, m_componentNames, m_componentMolarWeight,
                                       getLogLevel() > 0 && logger::internal::rank==0 );
}

template< typename PHASE >
void ReactiveBrineFluid< PHASE >::checkTablesParameters( real64 const pressure,
                                                         real64 const temperature ) const
{
  if( !m_checkPVTTablesRanges )
  {
    return;
  }

  real64 const temperatureInCelsius = units::convertKToC( temperature );
  try
  {
    m_phase->density.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase->viscosity.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase->enthalpy.checkTablesParameters( pressure, temperatureInCelsius );
  } catch( SimulationError const & ex )
  {
    string const errorMsg = GEOS_FMT( "Table input error (in table from {}).\n",
                                      stringutilities::join( m_phasePVTParaFiles ) );
    throw SimulationError( ex, errorMsg );
  }
}

template< typename PHASE >
typename ReactiveBrineFluid< PHASE > ::KernelWrapper
ReactiveBrineFluid< PHASE > ::createKernelWrapper()
{
  return KernelWrapper( *m_phase,
                        m_componentMolarWeight.toViewConst(),
                        m_useMass,
                        isThermal(),
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView(),
                        m_numPrimarySpecies,
                        *m_equilibriumReactions,
                        *m_kineticReactions,
                        m_primarySpeciesConcentration.toView(),
                        m_secondarySpeciesConcentration.toView(),
                        m_primarySpeciesTotalConcentration.toView(),
                        m_kineticReactionRates.toView() );
}

template< typename PHASE >
ReactiveBrineFluid< PHASE > ::KernelWrapper::
  KernelWrapper( PHASE const & phase,
                 arrayView1d< real64 const > componentMolarWeight,
                 bool const useMass,
                 bool const isThermal,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseProp::ViewType phaseEnthalpy,
                 PhaseProp::ViewType phaseInternalEnergy,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity,
                 integer const numPrimarySpecies,
                 chemicalReactions::EquilibriumReactions const & equilibriumReactions,
                 chemicalReactions::KineticReactions const & kineticReactions,
                 arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesConcentration,
                 arrayView2d< real64, compflow::USD_COMP > const & secondarySpeciesConcentration,
                 arrayView2d< real64, compflow::USD_COMP > const & primarySpeciesTotalConcentration,
                 arrayView2d< real64, compflow::USD_COMP > const & kineticReactionRates )
  : ReactiveMultiFluid::KernelWrapper( std::move( componentMolarWeight ),
                                       useMass,
                                       std::move( phaseFraction ),
                                       std::move( phaseDensity ),
                                       std::move( phaseMassDensity ),
                                       std::move( phaseViscosity ),
                                       std::move( phaseEnthalpy ),
                                       std::move( phaseInternalEnergy ),
                                       std::move( phaseCompFraction ),
                                       std::move( totalDensity ),
                                       numPrimarySpecies,
                                       equilibriumReactions,
                                       kineticReactions,
                                       primarySpeciesConcentration,
                                       secondarySpeciesConcentration,
                                       primarySpeciesTotalConcentration,
                                       kineticReactionRates ),
  m_isThermal( isThermal ),
  m_phase( phase.createKernelWrapper() )
{}

// explicit instantiation of the model template; unfortunately we can't use the aliases for this
template class ReactiveBrineFluid< PhaseModel< PVTProps::WaterDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction > >;
template class ReactiveBrineFluid< PhaseModel< PVTProps::WaterDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy > >;


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveBrine, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveBrineThermal, string const &, Group * const )

} //namespace constitutive

} //namespace geos
