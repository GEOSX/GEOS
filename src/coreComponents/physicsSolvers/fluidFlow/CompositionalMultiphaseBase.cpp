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
 * @file CompositionalMultiphaseBase.cpp
 */

#include "CompositionalMultiphaseBase.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/diffusion/DiffusionFields.hpp"
#include "constitutive/diffusion/DiffusionSelector.hpp"
#include "constitutive/dispersion/DispersionFields.hpp"
#include "constitutive/dispersion/DispersionSelector.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"
#include "constitutive/solid/SolidInternalEnergy.hpp"
#include "constitutive/thermalConductivity/MultiPhaseThermalConductivitySelector.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "physicsSolvers/fluidFlow/SourceFluxStatistics.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/IsothermalCompositionalMultiphaseFVMKernels.hpp" // TODO this should not be here
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalCompositionalMultiphaseBaseKernels.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

CompositionalMultiphaseBase::CompositionalMultiphaseBase( const string & name,
                                                          Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_hasCapPressure( 0 ),
  m_hasDiffusion( 0 ),
  m_hasDispersion( 0 ),
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 ),
  m_useTotalMassEquation( 1 ),
  m_useSimpleAccumulation( 1 ),
  m_minCompDens( isothermalCompositionalMultiphaseBaseKernels::minDensForDivision )
{
//START_SPHINX_INCLUDE_00
  this->registerWrapper( viewKeyStruct::inputTemperatureString(), &m_inputTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Temperature" );
//END_SPHINX_INCLUDE_00

  // TODO: add more description on how useMass can alter the simulation (convergence issues?). How does it interact with wells useMass?
  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( GEOS_FMT( "Use mass formulation instead of molar. Warning : Affects {} rates units.",
                              SourceFluxBoundaryCondition::catalogName() ) );

  this->registerWrapper( viewKeyStruct::solutionChangeScalingFactorString(), &m_solutionChangeScalingFactor ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5 ).
    setDescription( "Damping factor for solution change targets" );
  this->registerWrapper( viewKeyStruct::targetRelativePresChangeString(), &m_targetRelativePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.2 ).
    setDescription( "Target (relative) change in pressure in a time step (expected value between 0 and 1)" );
  this->registerWrapper( viewKeyStruct::targetRelativeTempChangeString(), &m_targetRelativeTempChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.2 ).
    setDescription( "Target (relative) change in temperature in a time step (expected value between 0 and 1)" );
  this->registerWrapper( viewKeyStruct::targetPhaseVolFracChangeString(), &m_targetPhaseVolFracChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.2 ).
    setDescription( "Target (absolute) change in phase volume fraction in a time step" );
  this->registerWrapper( viewKeyStruct::targetRelativeCompDensChangeString(), &m_targetRelativeCompDensChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( LvArray::NumericLimits< real64 >::max ). // disabled by default
    setDescription( "Target (relative) change in component density in a time step" );

  this->registerWrapper( viewKeyStruct::maxCompFracChangeString(), &m_maxCompFracChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5 ).
    setDescription( "Maximum (absolute) change in a component fraction in a Newton iteration" );
  this->registerWrapper( viewKeyStruct::maxRelativePresChangeString(), &m_maxRelativePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5 ).
    setDescription( "Maximum (relative) change in pressure in a Newton iteration" );
  this->registerWrapper( viewKeyStruct::maxRelativeTempChangeString(), &m_maxRelativeTempChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5 ).
    setDescription( "Maximum (relative) change in temperature in a Newton iteration" );
  this->registerWrapper( viewKeyStruct::maxRelativeCompDensChangeString(), &m_maxRelativeCompDensChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( LvArray::NumericLimits< real64 >::max/1.0e100 ). // disabled by default
    setDescription( "Maximum (relative) change in a component density in a Newton iteration" );

  this->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString(), &m_allowCompDensChopping ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative compositions is allowed" );

  this->registerWrapper( viewKeyStruct::targetFlowCFLString(), &m_targetFlowCFL ).
    setApplyDefaultValue( -1. ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Target CFL condition `CFL condition <http://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition>`_"
                    " when computing the next timestep." );

  this->registerWrapper( viewKeyStruct::useTotalMassEquationString(), &m_useTotalMassEquation ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether total mass equation is used" );

  this->registerWrapper( viewKeyStruct::useSimpleAccumulationString(), &m_useSimpleAccumulation ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether simple accumulation form is used" );

  this->registerWrapper( viewKeyStruct::minCompDensString(), &m_minCompDens ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( isothermalCompositionalMultiphaseBaseKernels::minDensForDivision ).
    setDescription( "Minimum allowed global component density" );

  this->registerWrapper( viewKeyStruct::maxSequentialCompDensChangeString(), &m_maxSequentialCompDensChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (absolute) component density change in a sequential iteration, used for outer loop convergence check" );

  this->registerWrapper( viewKeyStruct::minScalingFactorString(), &m_minScalingFactor ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.01 ).
    setDescription( "Minimum value for solution scaling factor" );
}

void CompositionalMultiphaseBase::postInputInitialization()
{
  FlowSolverBase::postInputInitialization();

  GEOS_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                        getWrapperDataContext( viewKeyStruct::maxCompFracChangeString() ) <<
                        ": The maximum absolute change in component fraction in a Newton iteration must be smaller or equal to 1.0" );
  GEOS_ERROR_IF_LE_MSG( m_maxCompFracChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::maxCompFracChangeString() ) <<
                        ": The maximum absolute change in component fraction in a Newton iteration must be larger than 0.0" );
  GEOS_ERROR_IF_LE_MSG( m_maxRelativePresChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::maxRelativePresChangeString() ) <<
                        ": The maximum relative change in pressure in a Newton iteration must be larger than 0.0" );
  GEOS_ERROR_IF_LE_MSG( m_maxRelativeTempChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::maxRelativeTempChangeString() ) <<
                        ": The maximum relative change in temperature in a Newton iteration must be larger than 0.0" );
  GEOS_ERROR_IF_LE_MSG( m_maxRelativeCompDensChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::maxRelativeCompDensChangeString() ) <<
                        ": The maximum relative change in component density in a Newton iteration must be larger than 0.0" );
  GEOS_ERROR_IF_LE_MSG( m_targetRelativePresChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::targetRelativePresChangeString() ) <<
                        ": The target relative change in pressure in a time step must be larger than 0.0" );
  GEOS_ERROR_IF_LE_MSG( m_targetRelativeTempChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::targetRelativeTempChangeString() ) <<
                        ": The target relative change in temperature in a time step must be larger than 0.0" );
  GEOS_ERROR_IF_LE_MSG( m_targetPhaseVolFracChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::targetPhaseVolFracChangeString() ) <<
                        ": The target change in phase volume fraction in a time step must be larger than 0.0" );
  GEOS_ERROR_IF_LE_MSG( m_targetRelativeCompDensChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::targetPhaseVolFracChangeString() ) <<
                        ": The target change in component density in a time step must be larger than 0.0" );

  GEOS_ERROR_IF_LT_MSG( m_solutionChangeScalingFactor, 0.0,
                        getWrapperDataContext( viewKeyStruct::solutionChangeScalingFactorString() ) <<
                        ": The solution change scaling factor must be larger or equal to 0.0" );
  GEOS_ERROR_IF_GT_MSG( m_solutionChangeScalingFactor, 1.0,
                        getWrapperDataContext( viewKeyStruct::solutionChangeScalingFactorString() ) <<
                        ": The solution change scaling factor must be smaller or equal to 1.0" );
  GEOS_ERROR_IF_LE_MSG( m_minScalingFactor, 0.0,
                        getWrapperDataContext( viewKeyStruct::minScalingFactorString() ) <<
                        ": The minumum scaling factor must be larger than 0.0" );
  GEOS_ERROR_IF_GT_MSG( m_minScalingFactor, 1.0,
                        getWrapperDataContext( viewKeyStruct::minScalingFactorString() ) <<
                        ": The minumum scaling factor must be smaller or equal to 1.0" );

  if( m_isThermal && m_useSimpleAccumulation == 1 ) // useSimpleAccumulation is not yet compatible with thermal
  {
    GEOS_LOG_RANK_0( "'useSimpleAccumulation' is not yet implemented for thermal simulation. Switched to phase sum accumulation." );
    m_useSimpleAccumulation = 0;
  }
}

void CompositionalMultiphaseBase::registerDataOnMesh( Group & meshBodies )
{
  using namespace fields::flow;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 0. Find a "reference" fluid model name (at this point, models are already attached to subregions)
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      if( m_referenceFluidModelName.empty() )
      {
        m_referenceFluidModelName = getConstitutiveName< MultiFluidBase >( subRegion );
      }

      // If at least one region has a capillary pressure model, consider it enabled for all
      string const capPresName = getConstitutiveName< CapillaryPressureBase >( subRegion );
      if( !capPresName.empty() )
      {
        m_hasCapPressure = true;
      }

      // If at least one region has a diffusion model, consider it enabled for all
      string const diffusionName = getConstitutiveName< DiffusionBase >( subRegion );
      if( !diffusionName.empty() )
      {
        m_hasDiffusion = true;
      }

      // If at least one region has a dispersion model, consider it enabled for all
      string const dispersionName = getConstitutiveName< DispersionBase >( subRegion );
      if( !dispersionName.empty() )
      {
        GEOS_ERROR( "Dispersion is not supported yet, please remove this model from this XML file" );
        m_hasDispersion = true;
      }

    } );
  } );

  // 1. Set key dimensions of the problem
  // Check needed to avoid errors when running in schema generation mode.
  if( !m_referenceFluidModelName.empty() )
  {
    MultiFluidBase const & referenceFluid = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );
    m_numPhases = referenceFluid.numFluidPhases();
    m_numComponents = referenceFluid.numFluidComponents();
    m_isThermal = referenceFluid.isThermal();
  }

  // n_c components + one pressure ( + one temperature if needed )
  m_numDofPerCell = m_isThermal ? m_numComponents + 2 : m_numComponents + 1;

  // 2. Register and resize all fields as necessary
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      if( m_hasCapPressure )
      {
        subRegion.registerWrapper< string >( viewKeyStruct::capPressureNamesString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRestartFlags( RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 ).
          setDescription( "Name of the capillary pressure constitutive model to use" ).
          reference();

        string & capPresName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        capPresName = getConstitutiveName< CapillaryPressureBase >( subRegion );
        GEOS_THROW_IF( capPresName.empty(),
                       GEOS_FMT( "{}: Capillary pressure model not found on subregion {}",
                                 getDataContext(), subRegion.getDataContext() ),
                       InputError );
      }

      if( m_hasDiffusion )
      {
        subRegion.registerWrapper< string >( viewKeyStruct::diffusionNamesString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRestartFlags( RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 ).
          setDescription( "Name of the diffusion constitutive model to use" );

        string & diffusionName = subRegion.getReference< string >( viewKeyStruct::diffusionNamesString() );
        diffusionName = getConstitutiveName< DiffusionBase >( subRegion );
        GEOS_THROW_IF( diffusionName.empty(),
                       GEOS_FMT( "Diffusion model not found on subregion {}", subRegion.getName() ),
                       InputError );
      }

      if( m_hasDispersion )
      {
        subRegion.registerWrapper< string >( viewKeyStruct::dispersionNamesString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRestartFlags( RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 ).
          setDescription( "Name of the dispersion constitutive model to use" );

        string & dispersionName = subRegion.getReference< string >( viewKeyStruct::dispersionNamesString() );
        dispersionName = getConstitutiveName< DispersionBase >( subRegion );
        GEOS_THROW_IF( dispersionName.empty(),
                       GEOS_FMT( "Dispersion model not found on subregion {}", subRegion.getName() ),
                       InputError );
      }

      if( m_targetFlowCFL > 0 )
      {
        subRegion.registerField< fields::flow::phaseOutflux >( getName() ).
          reference().resizeDimension< 1 >( m_numPhases );
        subRegion.registerField< fields::flow::componentOutflux >( getName() ).
          reference().resizeDimension< 1 >( m_numComponents );
        subRegion.registerField< fields::flow::phaseCFLNumber >( getName() );
        subRegion.registerField< fields::flow::componentCFLNumber >( getName() );
      }

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      subRegion.registerField< pressureScalingFactor >( getName() );
      subRegion.registerField< temperatureScalingFactor >( getName() );
      subRegion.registerField< globalCompDensityScalingFactor >( getName() );

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.

      subRegion.registerField< globalCompDensity >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< globalCompDensity_n >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );
      if( m_isFixedStressPoromechanicsUpdate )
      {
        subRegion.registerField< globalCompDensity_k >( getName() ).
          setDimLabels( 1, fluid.componentNames() ).
          reference().resizeDimension< 1 >( m_numComponents );
      }

      subRegion.registerField< globalCompFraction >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< dGlobalCompFraction_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numComponents, m_numComponents );

      subRegion.registerField< phaseVolumeFraction >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< dPhaseVolumeFraction >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents + 2 ); // dP, dT, dC

      subRegion.registerField< phaseMobility >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< dPhaseMobility >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents + 2 ); // dP, dT, dC

      // needed for time step selector
      subRegion.registerField< phaseVolumeFraction_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< compAmount >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< compAmount_n >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );

    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      // We make the assumption that component names are uniform across the fluid models used in the simulation
      MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );

      // TODO: add conditional registration later, this is only needed when there is a face-based Dirichlet BC
      faceManager.registerField< facePressure >( getName() );
      faceManager.registerField< faceTemperature >( getName() );
      faceManager.registerField< faceGlobalCompFraction >( getName() ).
        setDimLabels( 1, fluid0.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
    }
  } );
}

void CompositionalMultiphaseBase::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
  GEOS_THROW_IF( fluidName.empty(),
                 GEOS_FMT( "{}: Fluid model not found on subregion {}",
                           getDataContext(), subRegion.getDataContext() ),
                 InputError );

  string & relPermName = subRegion.registerWrapper< string >( viewKeyStruct::relPermNamesString() ).
                           setPlotLevel( PlotLevel::NOPLOT ).
                           setRestartFlags( RestartFlags::NO_WRITE ).
                           setSizedFromParent( 0 ).
                           setDescription( "Name of the relative permeability constitutive model to use" ).
                           reference();

  relPermName = getConstitutiveName< RelativePermeabilityBase >( subRegion );

  GEOS_THROW_IF( relPermName.empty(),
                 GEOS_FMT( "{}: Relative permeability model not found on subregion {}",
                           getDataContext(), subRegion.getDataContext() ),
                 InputError );

  if( m_isThermal )
  {
    string & thermalConductivityName = subRegion.registerWrapper< string >( viewKeyStruct::thermalConductivityNamesString() ).
                                         setPlotLevel( PlotLevel::NOPLOT ).
                                         setRestartFlags( RestartFlags::NO_WRITE ).
                                         setSizedFromParent( 0 ).
                                         setDescription( "Name of the thermal conductivity constitutive model to use" ).
                                         reference();

    thermalConductivityName = getConstitutiveName< MultiPhaseThermalConductivityBase >( subRegion );
    GEOS_THROW_IF( thermalConductivityName.empty(),
                   GEOS_FMT( "{}: Thermal conductivity model not found on subregion {}",
                             getDataContext(), subRegion.getDataContext() ),
                   InputError );
  }
}


namespace
{

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOS_THROW_IF_NE_MSG( lhs.numFluidPhases(), rhs.numFluidPhases(),
                        GEOS_FMT( "Mismatch in number of phases between constitutive models {} and {}",
                                  lhs.getDataContext(), rhs.getDataContext() ),
                        InputError );

  for( integer ip = 0; ip < lhs.numFluidPhases(); ++ip )
  {
    GEOS_THROW_IF_NE_MSG( lhs.phaseNames()[ip], rhs.phaseNames()[ip],
                          GEOS_FMT( "Mismatch in phase names between constitutive models {} and {}",
                                    lhs.getDataContext(), rhs.getDataContext() ),
                          InputError );
  }
}

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOS_THROW_IF_NE_MSG( lhs.numFluidComponents(), rhs.numFluidComponents(),
                        GEOS_FMT( "Mismatch in number of components between constitutive models {} and {}",
                                  lhs.getDataContext(), rhs.getDataContext() ),
                        InputError );

  for( integer ic = 0; ic < lhs.numFluidComponents(); ++ic )
  {
    GEOS_THROW_IF_NE_MSG( lhs.componentNames()[ic], rhs.componentNames()[ic],
                          GEOS_FMT( "Mismatch in component names between constitutive models {} and {}",
                                    lhs.getDataContext(), rhs.getDataContext() ),
                          InputError );
  }
}

}

void CompositionalMultiphaseBase::initializeAquiferBC( ConstitutiveManager const & cm ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition & bc )
  {
    MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );

    // set the gravity vector (needed later for the potential diff calculations)
    bc.setGravityVector( gravityVector() );

    // set the water phase index in the Aquifer boundary condition
    // note: if the water phase is not found, the fluid model is going to throw an error
    integer const waterPhaseIndex = fluid0.getWaterPhaseIndex();
    bc.setWaterPhaseIndex( waterPhaseIndex );

    arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac = bc.getWaterPhaseComponentFraction();
    arrayView1d< string const > const & aquiferWaterPhaseCompNames = bc.getWaterPhaseComponentNames();

    GEOS_ERROR_IF_NE_MSG( fluid0.numFluidComponents(), aquiferWaterPhaseCompFrac.size(),
                          getDataContext() << ": Mismatch in number of components between constitutive model "
                                           << fluid0.getName() << " and the water phase composition in aquifer " << bc.getName() );

    for( integer ic = 0; ic < fluid0.numFluidComponents(); ++ic )
    {
      GEOS_ERROR_IF_NE_MSG( fluid0.componentNames()[ic], aquiferWaterPhaseCompNames[ic],
                            getDataContext() << ": Mismatch in component names between constitutive model "
                                             << fluid0.getName() << " and the water phase components in aquifer " << bc.getName() );
    }
  } );
}


void CompositionalMultiphaseBase::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 1. Validate various models against each other (must have same phases and components)
  validateConstitutiveModels( domain );

  // 2. Set the value of temperature
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const temp = subRegion.getField< fields::flow::temperature >();
      temp.setValues< parallelHostPolicy >( m_inputTemperature );
    } );
  } );

  // 3. Initialize and validate the aquifer boundary condition
  initializeAquiferBC( cm );

}

void CompositionalMultiphaseBase::validateConstitutiveModels( DomainPartition const & domain ) const
{
  GEOS_MARK_FUNCTION;

  ConstitutiveManager const & cm = domain.getConstitutiveManager();
  MultiFluidBase const & referenceFluid = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      compareMultiphaseModels( fluid, referenceFluid );
      compareMulticomponentModels( fluid, referenceFluid );

      bool const isFluidModelThermal = fluid.isThermal();
      GEOS_THROW_IF( m_isThermal && !isFluidModelThermal,
                     GEOS_FMT( "CompositionalMultiphaseBase {}: the thermal option is enabled in the solver, but the fluid model {} is incompatible with the thermal option",
                               getDataContext(), fluid.getDataContext() ),
                     InputError );
      GEOS_THROW_IF( !m_isThermal && isFluidModelThermal,
                     GEOS_FMT( "CompositionalMultiphaseBase {}: the thermal option is enabled in fluid model {}, but the solver options are incompatible with the thermal option",
                               getDataContext(), fluid.getDataContext() ),
                     InputError );

      string const & relpermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relPerm = getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
      compareMultiphaseModels( relPerm, referenceFluid );

      if( m_hasCapPressure )
      {
        string const & capPressureName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        CapillaryPressureBase const & capPressure = getConstitutiveModel< CapillaryPressureBase >( subRegion, capPressureName );
        compareMultiphaseModels( capPressure, referenceFluid );
      }

      if( m_hasDiffusion )
      {
        string const & diffusionName = subRegion.getReference< string >( viewKeyStruct::diffusionNamesString() );
        DiffusionBase const & diffusion = getConstitutiveModel< DiffusionBase >( subRegion, diffusionName );
        compareMultiphaseModels( diffusion, referenceFluid );
      }

      if( m_isThermal )
      {
        string const & thermalConductivityName = subRegion.getReference< string >( viewKeyStruct::thermalConductivityNamesString() );
        MultiPhaseThermalConductivityBase const & conductivity = getConstitutiveModel< MultiPhaseThermalConductivityBase >( subRegion, thermalConductivityName );
        compareMultiphaseModels( conductivity, referenceFluid );
      }
    } );
  } );
}

void CompositionalMultiphaseBase::updateGlobalComponentFraction( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  isothermalCompositionalMultiphaseBaseKernels::
    GlobalComponentFractionKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               dataGroup );

}

real64 CompositionalMultiphaseBase::updatePhaseVolumeFraction( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, fluidName );

  real64 maxDeltaPhaseVolFrac  =
    m_isThermal ?
    thermalCompositionalMultiphaseBaseKernels::
      PhaseVolumeFractionKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 dataGroup,
                                                 fluid )
:    isothermalCompositionalMultiphaseBaseKernels::
      PhaseVolumeFractionKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 dataGroup,
                                                 fluid );

  return maxDeltaPhaseVolFrac;
}

void CompositionalMultiphaseBase::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const temp = dataGroup.getField< fields::flow::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
    dataGroup.getField< fields::flow::globalCompFraction >();

  string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, fluidName );

  constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    thermalCompositionalMultiphaseBaseKernels::
      FluidUpdateKernel::
      launch< ExecPolicy >( dataGroup.size(),
                            fluidWrapper,
                            pres,
                            temp,
                            compFrac );
  } );
}

void CompositionalMultiphaseBase::updateRelPermModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
    dataGroup.getField< fields::flow::phaseVolumeFraction >();

  string const & relPermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase & relPerm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relPermName );

  constitutive::constitutiveUpdatePassThru( relPerm, [&] ( auto & castedRelPerm )
  {
    typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();

    isothermalCompositionalMultiphaseBaseKernels::
      RelativePermeabilityUpdateKernel::
      launch< parallelDevicePolicy<> >( dataGroup.size(),
                                        relPermWrapper,
                                        phaseVolFrac );
  } );
}

void CompositionalMultiphaseBase::updateCapPressureModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  if( m_hasCapPressure )
  {
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      dataGroup.getField< fields::flow::phaseVolumeFraction >();

    string const & cappresName = dataGroup.getReference< string >( viewKeyStruct::capPressureNamesString() );
    CapillaryPressureBase & capPressure = getConstitutiveModel< CapillaryPressureBase >( dataGroup, cappresName );

    constitutive::constitutiveUpdatePassThru( capPressure, [&] ( auto & castedCapPres )
    {
      typename TYPEOFREF( castedCapPres ) ::KernelWrapper capPresWrapper = castedCapPres.createKernelWrapper();

      isothermalCompositionalMultiphaseBaseKernels::
        CapillaryPressureUpdateKernel::
        launch< parallelDevicePolicy<> >( dataGroup.size(),
                                          capPresWrapper,
                                          phaseVolFrac );
    } );
  }
}

void CompositionalMultiphaseBase::updateCompAmount( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  string const & solidName = subRegion.template getReference< string >( viewKeyStruct::solidNamesString() );
  CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
  arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();
  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::flow::globalCompDensity >();
  arrayView2d< real64, compflow::USD_COMP > const compAmount = subRegion.getField< fields::flow::compAmount >();

  integer const numComp = m_numComponents;

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( integer ic = 0; ic < numComp; ++ic )
    {
      compAmount[ei][ic] = porosity[ei][0] * volume[ei] * compDens[ei][ic];
    }
  } );
}

void CompositionalMultiphaseBase::updateEnergy( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  string const & solidName = subRegion.template getReference< string >( viewKeyStruct::solidNamesString() );
  CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
  arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();
  arrayView2d< real64 const > rockInternalEnergy = porousMaterial.getInternalEnergy();
  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac = subRegion.getField< fields::flow::phaseVolumeFraction >();
  string const & fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
  arrayView3d< real64 const, multifluid::USD_PHASE > const phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const, multifluid::USD_PHASE > const phaseInternalEnergy = fluid.phaseInternalEnergy();

  arrayView1d< real64 > const energy = subRegion.getField< fields::flow::energy >();

  integer const numPhases = m_numPhases;

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    energy[ei] = volume[ei] * (1.0 - porosity[ei][0]) * rockInternalEnergy[ei][0];
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      energy[ei] += volume[ei] * porosity[ei][0] * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip] * phaseInternalEnergy[ei][0][ip];
    }
  } );
}

void CompositionalMultiphaseBase::updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const
{
  arrayView1d< real64 const > const temp = dataGroup.getField< fields::flow::temperature >();

  string const & solidInternalEnergyName = dataGroup.getReference< string >( viewKeyStruct::solidInternalEnergyNamesString() );
  SolidInternalEnergy & solidInternalEnergy = getConstitutiveModel< SolidInternalEnergy >( dataGroup, solidInternalEnergyName );

  SolidInternalEnergy::KernelWrapper solidInternalEnergyWrapper = solidInternalEnergy.createKernelUpdates();

  // TODO: this should go somewhere, handle the case of flow in fracture, etc

  thermalCompositionalMultiphaseBaseKernels::
    SolidInternalEnergyUpdateKernel::
    launch< parallelDevicePolicy<> >( dataGroup.size(),
                                      solidInternalEnergyWrapper,
                                      temp );
}

real64 CompositionalMultiphaseBase::updateFluidState( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  updateGlobalComponentFraction( subRegion );
  updateFluidModel( subRegion );
  updateCompAmount( subRegion );
  real64 const maxDeltaPhaseVolFrac = updatePhaseVolumeFraction( subRegion );
  updateRelPermModel( subRegion );
  updatePhaseMobility( subRegion );
  updateCapPressureModel( subRegion );

  // note1: for now, thermal conductivity is treated explicitly, so no update here
  // note2: for now, diffusion and dispersion are also treated explicitly
  return maxDeltaPhaseVolFrac;
}

void CompositionalMultiphaseBase::initializeFluidState( MeshLevel & mesh,
                                                        DomainPartition & domain,
                                                        arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;

  integer const numComp = m_numComponents;

  // 1. Compute hydrostatic equilibrium in the regions for which corresponding field specification tag has been specified
  computeHydrostaticEquilibrium();

  mesh.getElemManager().forElementSubRegions( regionNames,
                                              [&]( localIndex const,
                                                   ElementSubRegionBase & subRegion )
  {
    // 2. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    updateFluidModel( subRegion );

    // 3. Back-calculate global component densities from fractions and total fluid density
    // in order to initialize the primary solution variables
    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
    arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();

    arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
      subRegion.getField< fields::flow::globalCompFraction >();
    arrayView2d< real64, compflow::USD_COMP > const compDens =
      subRegion.getField< fields::flow::globalCompDensity >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      for( integer ic = 0; ic < numComp; ++ic )
      {
        compDens[ei][ic] = totalDens[ei][0] * compFrac[ei][ic];
      }
    } );
  } );

  // with initial component densities defined - check if they need to be corrected to avoid zero diags etc
  chopNegativeDensities( domain );

  // for some reason CUDA does not want the host_device lambda to be defined inside the generic lambda
  // I need the exact type of the subRegion for updateSolidflowProperties to work well.
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                              SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                           auto & subRegion )
  {
    // 4. Initialize/update dependent state quantities

    // 4.1 Update the constitutive models that only depend on
    //      - the primary variables
    //      - the fluid constitutive quantities (as they have already been updated)
    // We postpone the other constitutive models for now
    // In addition, to avoid multiplying permeability/porosity bay netToGross in the assembly kernel, we do it once and for all here
    arrayView1d< real64 const > const netToGross = subRegion.template getField< fields::flow::netToGross >();
    CoupledSolidBase const & porousSolid =
      getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
    PermeabilityBase const & permeabilityModel =
      getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() ) );
    permeabilityModel.scaleHorizontalPermeability( netToGross );
    porousSolid.scaleReferencePorosity( netToGross );
    saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
    updatePorosityAndPermeability( subRegion );
    updateCompAmount( subRegion );
    updatePhaseVolumeFraction( subRegion );

    // Now, we initialize and update each constitutive model one by one

    // 4.2 Save the computed porosity into the old porosity
    //
    // Note:
    // - This must be called after updatePorosityAndPermeability
    // - This step depends on porosity
    string const & solidName = subRegion.template getReference< string >( viewKeyStruct::solidNamesString() );
    CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
    porousMaterial.initializeState();

    // 4.3 Initialize/update the relative permeability model using the initial phase volume fraction
    //     This is needed to handle relative permeability hysteresis
    //     Also, initialize the fluid model
    //
    // Note:
    // - This must be called after updatePhaseVolumeFraction
    // - This step depends on phaseVolFraction

    // initialized phase volume fraction
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      subRegion.template getField< fields::flow::phaseVolumeFraction >();

    string const & relpermName = subRegion.template getReference< string >( viewKeyStruct::relPermNamesString() );
    RelativePermeabilityBase & relPermMaterial =
      getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
    relPermMaterial.saveConvergedPhaseVolFractionState( phaseVolFrac ); // this needs to happen before calling updateRelPermModel
    updateRelPermModel( subRegion );
    relPermMaterial.saveConvergedState(); // this needs to happen after calling updateRelPermModel

    string const & fluidName = subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() );
    MultiFluidBase & fluidMaterial = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
    fluidMaterial.initializeState();

    // 4.4 Then, we initialize/update the capillary pressure model
    //
    // Note:
    // - This must be called after updatePorosityAndPermeability
    // - This step depends on porosity and permeability
    if( m_hasCapPressure )
    {
      // initialized porosity
      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

      string const & permName = subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeabilityMaterial =
        getConstitutiveModel< PermeabilityBase >( subRegion, permName );
      // initialized permeability
      arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

      string const & capPressureName = subRegion.template getReference< string >( viewKeyStruct::capPressureNamesString() );
      CapillaryPressureBase const & capPressureMaterial =
        getConstitutiveModel< CapillaryPressureBase >( subRegion, capPressureName );
      capPressureMaterial.initializeRockState( porosity, permeability ); // this needs to happen before calling updateCapPressureModel
      updateCapPressureModel( subRegion );
    }

    // 4.5 Update the phase mobility
    //
    // Note:
    // - This must be called after updateRelPermModel
    // - This step depends phaseRelPerm
    updatePhaseMobility( subRegion );

    // 4.6 We initialize the rock thermal quantities: conductivity and solid internal energy
    //
    // Note:
    // - This must be called after updatePorosityAndPermeability and updatePhaseVolumeFraction
    // - This step depends on porosity and phaseVolFraction
    if( m_isThermal )
    {
      // initialized porosity
      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

      string const & thermalConductivityName = subRegion.template getReference< string >( viewKeyStruct::thermalConductivityNamesString() );
      MultiPhaseThermalConductivityBase const & conductivityMaterial =
        getConstitutiveModel< MultiPhaseThermalConductivityBase >( subRegion, thermalConductivityName );
      conductivityMaterial.initializeRockFluidState( porosity, phaseVolFrac );
      // note that there is nothing to update here because thermal conductivity is explicit for now

      updateSolidInternalEnergyModel( subRegion );
      string const & solidInternalEnergyName = subRegion.template getReference< string >( viewKeyStruct::solidInternalEnergyNamesString() );
      SolidInternalEnergy const & solidInternalEnergyMaterial =
        getConstitutiveModel< SolidInternalEnergy >( subRegion, solidInternalEnergyName );
      solidInternalEnergyMaterial.saveConvergedState();

      updateEnergy( subRegion );
    }

    // Step 4.7: if the diffusion and/or dispersion is/are supported, initialize the two models
    if( m_hasDiffusion )
    {
      string const & diffusionName = subRegion.template getReference< string >( viewKeyStruct::diffusionNamesString() );
      DiffusionBase const & diffusionMaterial = getConstitutiveModel< DiffusionBase >( subRegion, diffusionName );
      arrayView1d< real64 const > const temperature = subRegion.template getField< fields::flow::temperature >();
      diffusionMaterial.initializeTemperatureState( temperature );
    }
    if( m_hasDispersion )
    {
      string const & dispersionName = subRegion.template getReference< string >( viewKeyStruct::dispersionNamesString() );
      DispersionBase const & dispersionMaterial = getConstitutiveModel< DispersionBase >( subRegion, dispersionName );
      GEOS_UNUSED_VAR( dispersionMaterial );
      // TODO: compute the phase velocities here
      //dispersionMaterial.saveConvergedVelocitySate( phaseVelovity );
    }

  } );

  // 5. Save initial pressure
  mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 > const initPres = subRegion.getField< fields::flow::initialPressure >();
    arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
    arrayView1d< real64 > const initTemp = subRegion.template getField< fields::flow::initialTemperature >();
    initPres.setValues< parallelDevicePolicy<> >( pres );
    initTemp.setValues< parallelDevicePolicy<> >( temp );
  } );
}

void CompositionalMultiphaseBase::computeHydrostaticEquilibrium()
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  integer const numComps = m_numComponents;
  integer const numPhases = m_numPhases;

  real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  // Step 1: count individual equilibriums (there may be multiple ones)

  std::map< string, localIndex > equilNameToEquilId;
  localIndex equilCounter = 0;

  fsManager.forSubGroups< EquilibriumInitialCondition >( [&] ( EquilibriumInitialCondition const & bc )
  {

    // collect all the equilibrium names to idx
    equilNameToEquilId[bc.getName()] = equilCounter;
    equilCounter++;

    // check that the gravity vector is aligned with the z-axis
    GEOS_THROW_IF( !isZero( gravVector[0] ) || !isZero( gravVector[1] ),
                   getCatalogName() << " " << getDataContext() <<
                   ": the gravity vector specified in this simulation (" << gravVector[0] << " " << gravVector[1] << " " << gravVector[2] <<
                   ") is not aligned with the z-axis. \n"
                   "This is incompatible with the " << bc.getCatalogName() << " " << bc.getDataContext() <<
                   "used in this simulation. To proceed, you can either: \n" <<
                   "   - Use a gravityVector aligned with the z-axis, such as (0.0,0.0,-9.81)\n" <<
                   "   - Remove the hydrostatic equilibrium initial condition from the XML file",
                   InputError );

  } );

  if( equilCounter == 0 )
  {
    return;
  }

  // Step 2: find the min elevation and the max elevation in the targetSets

  array1d< real64 > globalMaxElevation( equilNameToEquilId.size() );
  array1d< real64 > globalMinElevation( equilNameToEquilId.size() );
  findMinMaxElevationInEquilibriumTarget( domain,
                                          equilNameToEquilId,
                                          globalMaxElevation,
                                          globalMinElevation );

  // Step 3: for each equil, compute a fine table with hydrostatic pressure vs elevation if the region is a target region

  // first compute the region filter
  std::set< string > regionFilter;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    for( string const & regionName : regionNames )
    {
      regionFilter.insert( regionName );
    }

    fsManager.apply< ElementSubRegionBase,
                     EquilibriumInitialCondition >( 0.0,
                                                    mesh,
                                                    EquilibriumInitialCondition::catalogName(),
                                                    [&] ( EquilibriumInitialCondition const & fs,
                                                          string const &,
                                                          SortedArrayView< localIndex const > const & targetSet,
                                                          ElementSubRegionBase & subRegion,
                                                          string const & )
    {
      // Step 3.1: retrieve the data necessary to construct the pressure table in this subregion

      integer const maxNumEquilIterations = fs.getMaxNumEquilibrationIterations();
      real64 const equilTolerance = fs.getEquilibrationTolerance();
      real64 const datumElevation = fs.getDatumElevation();
      real64 const datumPressure = fs.getDatumPressure();
      string const initPhaseName = fs.getInitPhaseName(); // will go away when GOC/WOC are implemented

      localIndex const equilIndex = equilNameToEquilId.at( fs.getName() );
      real64 const minElevation = LvArray::math::min( globalMinElevation[equilIndex], datumElevation );
      real64 const maxElevation = LvArray::math::max( globalMaxElevation[equilIndex], datumElevation );
      real64 const elevationIncrement = LvArray::math::min( fs.getElevationIncrement(), maxElevation - minElevation );
      localIndex const numPointsInTable = ( elevationIncrement > 0 ) ? std::ceil( (maxElevation - minElevation) / elevationIncrement ) + 1 : 1;

      real64 const eps = 0.1 * (maxElevation - minElevation); // we add a small buffer to only log in the pathological cases
      GEOS_LOG_RANK_0_IF( ( (datumElevation > globalMaxElevation[equilIndex]+eps)  || (datumElevation < globalMinElevation[equilIndex]-eps) ),
                          getCatalogName() << " " << getDataContext() <<
                          ": By looking at the elevation of the cell centers in this model, GEOS found that " <<
                          "the min elevation is " << globalMinElevation[equilIndex] << " and the max elevation is " <<
                          globalMaxElevation[equilIndex] << "\nBut, a datum elevation of " << datumElevation <<
                          " was specified in the input file to equilibrate the model.\n " <<
                          "The simulation is going to proceed with this out-of-bound datum elevation," <<
                          " but the initial condition may be inaccurate." );

      array1d< array1d< real64 > > elevationValues;
      array1d< real64 > pressureValues;
      elevationValues.resize( 1 );
      elevationValues[0].resize( numPointsInTable );
      pressureValues.resize( numPointsInTable );

      // Step 3.2: retrieve the user-defined tables (temperature and comp fraction)

      FunctionManager & functionManager = FunctionManager::getInstance();

      array1d< TableFunction::KernelWrapper > compFracTableWrappers;
      arrayView1d< string const > compFracTableNames = fs.getComponentFractionVsElevationTableNames();
      for( integer ic = 0; ic < numComps; ++ic )
      {
        TableFunction const & compFracTable = functionManager.getGroup< TableFunction >( compFracTableNames[ic] );
        compFracTableWrappers.emplace_back( compFracTable.createKernelWrapper() );
      }

      string const tempTableName = fs.getTemperatureVsElevationTableName();
      TableFunction const & tempTable = functionManager.getGroup< TableFunction >( tempTableName );
      TableFunction::KernelWrapper tempTableWrapper = tempTable.createKernelWrapper();

      // Step 3.3: retrieve the fluid model to compute densities
      // we end up with the same issue as in applyDirichletBC: there is not a clean way to retrieve the fluid info

      Group const & region = subRegion.getParent().getParent();
      auto itRegionFilter = regionFilter.find( region.getName() );
      if( itRegionFilter == regionFilter.end() )
      {
        return; // the region is not in target, there is nothing to do
      }
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      arrayView1d< string const > componentNames = fs.getComponentNames();
      GEOS_THROW_IF( fluid.componentNames().size() != componentNames.size(),
                     "Mismatch in number of components between constitutive model "
                     << fluid.getDataContext() << " and the Equilibrium initial condition " << fs.getDataContext(),
                     InputError );
      for( integer ic = 0; ic < fluid.numFluidComponents(); ++ic )
      {
        GEOS_THROW_IF( fluid.componentNames()[ic] != componentNames[ic],
                       "Mismatch in component names between constitutive model "
                       << fluid.getDataContext() << " and the Equilibrium initial condition " << fs.getDataContext(),
                       InputError );
      }

      // Note: for now, we assume that the reservoir is in a single-phase state at initialization
      arrayView1d< string const > phaseNames = fluid.phaseNames();
      auto const itPhaseNames = std::find( std::begin( phaseNames ), std::end( phaseNames ), initPhaseName );
      GEOS_THROW_IF( itPhaseNames == std::end( phaseNames ),
                     getCatalogName() << " " << getDataContext() << ": phase name " <<
                     initPhaseName << " not found in the phases of " << fluid.getDataContext(),
                     InputError );
      integer const ipInit = std::distance( std::begin( phaseNames ), itPhaseNames );

      // Step 3.4: compute the hydrostatic pressure values

      constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
      {
        using FluidType = TYPEOFREF( castedFluid );
        typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

        // note: inside this kernel, serialPolicy is used, and elevation/pressure values don't go to the GPU
        isothermalCompositionalMultiphaseBaseKernels::
          HydrostaticPressureKernel::ReturnType const returnValue =
          isothermalCompositionalMultiphaseBaseKernels::
            HydrostaticPressureKernel::launch( numPointsInTable,
                                               numComps,
                                               numPhases,
                                               ipInit,
                                               maxNumEquilIterations,
                                               equilTolerance,
                                               gravVector,
                                               minElevation,
                                               elevationIncrement,
                                               datumElevation,
                                               datumPressure,
                                               fluidWrapper,
                                               compFracTableWrappers.toViewConst(),
                                               tempTableWrapper,
                                               elevationValues.toNestedView(),
                                               pressureValues.toView() );

        GEOS_THROW_IF( returnValue ==  isothermalCompositionalMultiphaseBaseKernels::HydrostaticPressureKernel::ReturnType::FAILED_TO_CONVERGE,
                       getCatalogName() << " " << getDataContext() <<
                       ": hydrostatic pressure initialization failed to converge in region " << region.getName() << "! \n" <<
                       "Try to loosen the equilibration tolerance, or increase the number of equilibration iterations. \n" <<
                       "If nothing works, something may be wrong in the fluid model, see <Constitutive> ",
                       std::runtime_error );

        GEOS_LOG_RANK_0_IF( returnValue == isothermalCompositionalMultiphaseBaseKernels::HydrostaticPressureKernel::ReturnType::DETECTED_MULTIPHASE_FLOW,
                            getCatalogName() << " " << getDataContext() <<
                            ": currently, GEOS assumes that there is only one mobile phase when computing the hydrostatic pressure. \n" <<
                            "We detected multiple phases using the provided datum pressure, temperature, and component fractions. \n" <<
                            "Please make sure that only one phase is mobile at the beginning of the simulation. \n" <<
                            "If this is not the case, the problem will not be at equilibrium when the simulation starts" );

      } );

      // Step 3.5: create hydrostatic pressure table

      string const tableName = fs.getName() + "_" + subRegion.getName() + "_" + phaseNames[ipInit] + "_table";
      TableFunction * const presTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
      presTable->setTableCoordinates( elevationValues, { units::Distance } );
      presTable->setTableValues( pressureValues, units::Pressure );
      presTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
      TableFunction::KernelWrapper presTableWrapper = presTable->createKernelWrapper();

      // Step 4: assign pressure, temperature, and component fraction as a function of elevation
      // TODO: this last step should probably be delayed to wait for the creation of FaceElements
      // TODO: this last step should be modified to account for GOC and WOC
      arrayView2d< real64 const > const elemCenter =
        subRegion.getReference< array2d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

      arrayView1d< real64 > const pres = subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );
      arrayView1d< real64 > const temp = subRegion.getReference< array1d< real64 > >( fields::flow::temperature::key() );
      arrayView2d< real64, compflow::USD_COMP > const compFrac =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::globalCompFraction::key() );
      arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappersViewConst =
        compFracTableWrappers.toViewConst();

      RAJA::ReduceMin< parallelDeviceReduce, real64 > minPressure( LvArray::NumericLimits< real64 >::max );

      forAll< parallelDevicePolicy<> >( targetSet.size(), [targetSet,
                                                           elemCenter,
                                                           presTableWrapper,
                                                           tempTableWrapper,
                                                           compFracTableWrappersViewConst,
                                                           numComps,
                                                           minPressure,
                                                           pres,
                                                           temp,
                                                           compFrac] GEOS_HOST_DEVICE ( localIndex const i )
      {
        localIndex const k = targetSet[i];
        real64 const elevation = elemCenter[k][2];

        pres[k] = presTableWrapper.compute( &elevation );
        minPressure.min( pres[k] );
        temp[k] = tempTableWrapper.compute( &elevation );
        for( integer ic = 0; ic < numComps; ++ic )
        {
          compFrac[k][ic] = compFracTableWrappersViewConst[ic].compute( &elevation );
        }
      } );

      GEOS_ERROR_IF( minPressure.get() < 0.0,
                     GEOS_FMT( "{}: A negative pressure of {} Pa was found during hydrostatic initialization in region/subRegion {}/{}",
                               getDataContext(), minPressure.get(), region.getName(), subRegion.getName() ) );
    } );
  } );
}

void CompositionalMultiphaseBase::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::flow::pressure::key(),
                                       fields::flow::globalCompDensity::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames,
                                                                                                 [&]( localIndex const,
                                                                                                      auto & subRegion )
    {
      // set mass fraction flag on fluid models
      string const & fluidName = subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      fluid.setMassFlag( m_useMass );

      saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
      updatePorosityAndPermeability( subRegion );

      CoupledSolidBase const & porousSolid =
        getConstitutiveModel< CoupledSolidBase >( subRegion,
                                                  subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
      porousSolid.initializeState();
    } );

    // Initialize primary variables from applied initial conditions
    initializeFluidState( mesh, domain, regionNames );

    mesh.getElemManager().forElementRegions< SurfaceElementRegion >( regionNames,
                                                                     [&]( localIndex const,
                                                                          SurfaceElementRegion & region )
    {
      region.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
      {
        subRegion.getWrapper< real64_array >( fields::flow::hydraulicAperture::key() ).
          setApplyDefaultValue( region.getDefaultAperture() );
      } );
    } );

  } );

  // report to the user if some pore volumes are very small
  // note: this function is here because: 1) porosity has been initialized and 2) NTG has been applied
  validatePoreVolumes( domain );
}

void
CompositionalMultiphaseBase::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                real64 const & GEOS_UNUSED_PARAM( dt ),
                                                DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      saveConvergedState( subRegion );

      // update porosity, permeability
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updateFluidState( subRegion );
      // for thermal simulations, update solid internal energy
      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
      }

      // after the update, save the new saturation
      arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.template getField< fields::flow::phaseVolumeFraction >();
      arrayView2d< real64, compflow::USD_PHASE > const phaseVolFrac_n =
        subRegion.template getField< fields::flow::phaseVolumeFraction_n >();
      phaseVolFrac_n.setValues< parallelDevicePolicy<> >( phaseVolFrac );

    } );
  } );
}

void CompositionalMultiphaseBase::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  assembleAccumulationAndVolumeBalanceTerms( domain,
                                             dofManager,
                                             localMatrix,
                                             localRhs );

  if( m_isJumpStabilized )
  {
    assembleStabilizedFluxTerms( dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );
  }
  else
  {
    assembleFluxTerms( dt,
                       domain,
                       dofManager,
                       localMatrix,
                       localRhs );
  }
}

void CompositionalMultiphaseBase::assembleAccumulationAndVolumeBalanceTerms( DomainPartition & domain,
                                                                             DofManager const & dofManager,
                                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                             arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );

      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      if( m_isThermal )
      {
        thermalCompositionalMultiphaseBaseKernels::
          ElementBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     m_useTotalMassEquation,
                                                     dofKey,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     localMatrix,
                                                     localRhs );
      }
      else
      {
        isothermalCompositionalMultiphaseBaseKernels::
          ElementBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     m_useTotalMassEquation,
                                                     m_useSimpleAccumulation,
                                                     dofKey,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     localMatrix,
                                                     localRhs );
      }
    } );
  } );
}

void CompositionalMultiphaseBase::applyBoundaryConditions( real64 const time_n,
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  if( m_keepVariablesConstantDuringInitStep )
  {
    // this function is going to force the current flow state to be constant during the time step
    // this is used when the poromechanics solver is performing the stress initialization
    // TODO: in the future, a dedicated poromechanics kernel should eliminate the flow vars to construct a reduced system
    //       which will remove the need for this brittle passing aroung of flag
    keepVariablesConstantDuringInitStep( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
  }
  else
  {
    // apply pressure boundary conditions.
    applyDirichletBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

    // apply flux boundary conditions
    applySourceFluxBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

    // apply aquifer boundary conditions
    applyAquiferBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
  }
}

namespace
{
char const bcLogMessage[] =
  "CompositionalMultiphaseBase {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target elements (including ghost elements) is {}. "
  "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set.";
}

void CompositionalMultiphaseBase::applySourceFluxBC( real64 const time,
                                                     real64 const dt,
                                                     DofManager const & dofManager,
                                                     DomainPartition & domain,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Step 1: count individual source flux boundary conditions

  std::map< string, localIndex > bcNameToBcId;
  localIndex bcCounter = 0;

  fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&] ( SourceFluxBoundaryCondition const & bc )
  {
    // collect all the bc names to idx
    bcNameToBcId[bc.getName()] = bcCounter;
    bcCounter++;
  } );

  if( bcCounter == 0 )
  {
    return;
  }

  // Step 2: count the set size for each source flux (each source flux may have multiple target sets)

  array1d< globalIndex > bcAllSetsSize( bcNameToBcId.size() );

  computeSourceFluxSizeScalingFactor( time,
                                      dt,
                                      domain,
                                      bcNameToBcId,
                                      bcAllSetsSize.toView() );

  // Step 3: we are ready to impose the boundary condition, normalized by the set size

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    fsManager.apply< ElementSubRegionBase,
                     SourceFluxBoundaryCondition >( time + dt,
                                                    mesh,
                                                    SourceFluxBoundaryCondition::catalogName(),
                                                    [&]( SourceFluxBoundaryCondition const & fs,
                                                         string const & setName,
                                                         SortedArrayView< localIndex const > const & targetSet,
                                                         ElementSubRegionBase & subRegion,
                                                         string const & )
    {
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOS_LOG_RANK_0( GEOS_FMT( bcLogMessage,
                                   getName(), time+dt, fs.getCatalogName(), fs.getName(),
                                   setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
      }

      if( targetSet.size() == 0 )
      {
        return;
      }
      if( !subRegion.hasWrapper( dofKey ) )
      {
        if( fs.getLogLevel() >= 1 )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: trying to apply SourceFlux, but its targetSet named '{}' intersects with non-simulated region named '{}'.",
                                   getDataContext(), setName, subRegion.getName() ) );
        }
        return;
      }

      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      // Step 3.1: get the values of the source boundary condition that need to be added to the rhs
      // We don't use FieldSpecificationBase::applyConditionToSystem here because we want to account for the row permutation used in the
      // compositional solvers

      array1d< globalIndex > dofArray( targetSet.size() );
      array1d< real64 > rhsContributionArray( targetSet.size() );
      arrayView1d< real64 > rhsContributionArrayView = rhsContributionArray.toView();
      localIndex const rankOffset = dofManager.rankOffset();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > massProd( 0.0 );

      // note that the dofArray will not be used after this step (simpler to use dofNumber instead)
      fs.computeRhsContribution< FieldSpecificationAdd,
                                 parallelDevicePolicy<> >( targetSet.toViewConst(),
                                                           time + dt,
                                                           dt,
                                                           subRegion,
                                                           dofNumber,
                                                           rankOffset,
                                                           localMatrix,
                                                           dofArray.toView(),
                                                           rhsContributionArrayView,
                                                           [] GEOS_HOST_DEVICE ( localIndex const )
      {
        return 0.0;
      } );

      // Step 3.2: we are ready to add the right-hand side contributions, taking into account our equation layout

      // get the normalizer
      real64 const sizeScalingFactor = bcAllSetsSize[bcNameToBcId.at( fs.getName())];

      integer const fluidComponentId = fs.getComponent();
      integer const numFluidComponents = m_numComponents;
      integer const useTotalMassEquation = m_useTotalMassEquation;
      forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                           targetSet,
                                                           rankOffset,
                                                           ghostRank,
                                                           fluidComponentId,
                                                           numFluidComponents,
                                                           useTotalMassEquation,
                                                           dofNumber,
                                                           rhsContributionArrayView,
                                                           localRhs,
                                                           massProd] GEOS_HOST_DEVICE ( localIndex const a )
      {
        // we need to filter out ghosts here, because targetSet may contain them
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        real64 const rhsValue = rhsContributionArrayView[a] / sizeScalingFactor; // scale the contribution by the sizeScalingFactor here!
        massProd += rhsValue;
        if( useTotalMassEquation > 0 )
        {
          // for all "fluid components", we add the value to the total mass balance equation
          globalIndex const totalMassBalanceRow = dofNumber[ei] - rankOffset;
          localRhs[totalMassBalanceRow] += rhsValue;
          if( fluidComponentId < numFluidComponents - 1 )
          {
            globalIndex const compMassBalanceRow = totalMassBalanceRow + fluidComponentId + 1; // component mass bal equations are shifted
            localRhs[compMassBalanceRow] += rhsValue;
          }
        }
        else
        {
          globalIndex const compMassBalanceRow = dofNumber[ei] - rankOffset + fluidComponentId;
          localRhs[compMassBalanceRow] += rhsValue;
        }
      } );

      SourceFluxStatsAggregator::forAllFluxStatWrappers( subRegion, fs.getName(),
                                                         [&]( SourceFluxStatsAggregator::WrappedStats & wrapper )
      {
        // set the new sub-region statistics for this timestep
        array1d< real64 > massProdArr{ m_numComponents };
        massProdArr[fluidComponentId] = massProd.get();
        wrapper.gatherTimeStepStats( time, dt, massProdArr.toViewConst(), targetSet.size() );
      } );
    } );
  } );
}

bool CompositionalMultiphaseBase::validateDirichletBC( DomainPartition & domain,
                                                       real64 const time ) const
{
  constexpr integer MAX_NC = MultiFluidBase::MAX_NUM_COMPONENTS;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  bool bcConsistent = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    // map: regionName -> subRegionName -> setName -> numComps to check pressure/comp are present consistent
    map< string, map< string, map< string, ComponentMask< MAX_NC > > > > bcPresCompStatusMap;
    // map: regionName -> subRegionName -> setName check to that temperature is present/consistent
    map< string, map< string, set< string > > > bcTempStatusMap;

    // 1. Check pressure Dirichlet BCs
    fsManager.apply< ElementSubRegionBase >( time,
                                             mesh,
                                             fields::flow::pressure::key(),
                                             [&]( FieldSpecificationBase const &,
                                                  string const & setName,
                                                  SortedArrayView< localIndex const > const &,
                                                  ElementSubRegionBase & subRegion,
                                                  string const & )
    {
      // Check whether pressure has already been applied to this set
      string const & subRegionName = subRegion.getName();
      string const & regionName = subRegion.getParent().getParent().getName();

      auto & subRegionSetMap = bcPresCompStatusMap[regionName][subRegionName];
      if( subRegionSetMap.count( setName ) > 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( BCMessage::pressureConflict( regionName, subRegionName, setName,
                                                   fields::flow::pressure::key() ) );
      }
      subRegionSetMap[setName].setNumComp( m_numComponents );
    } );

    // 2. Check temperature Dirichlet BCs
    if( m_isThermal )
    {
      fsManager.apply< ElementSubRegionBase >( time,
                                               mesh,
                                               fields::flow::temperature::key(),
                                               [&]( FieldSpecificationBase const &,
                                                    string const & setName,
                                                    SortedArrayView< localIndex const > const &,
                                                    ElementSubRegionBase & subRegion,
                                                    string const & )
      {
        // Check whether temperature has already been applied to this set
        string const & subRegionName = subRegion.getName();
        string const & regionName = subRegion.getParent().getParent().getName();

        auto & tempSubRegionSetMap = bcTempStatusMap[regionName][subRegionName];
        if( tempSubRegionSetMap.count( setName ) > 0 )
        {
          bcConsistent = false;
          GEOS_WARNING( BCMessage::temperatureConflict( regionName, subRegionName, setName,
                                                        fields::flow::temperature::key() ) );
        }
        tempSubRegionSetMap.insert( setName );
      } );
    }

    // 3. Check composition BC (global component fraction)
    fsManager.apply< ElementSubRegionBase >( time,
                                             mesh,
                                             fields::flow::globalCompFraction::key(),
                                             [&] ( FieldSpecificationBase const & fs,
                                                   string const & setName,
                                                   SortedArrayView< localIndex const > const &,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      // 3.1 Check pressure, temperature, and record composition bc application
      string const & subRegionName = subRegion.getName(   );
      string const & regionName = subRegion.getParent().getParent().getName();
      integer const comp = fs.getComponent();

      auto & subRegionSetMap = bcPresCompStatusMap[regionName][subRegionName];
      if( subRegionSetMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( BCMessage::missingPressure( regionName, subRegionName, setName,
                                                  fields::flow::pressure::key() ) );
      }
      if( m_isThermal )
      {
        auto & tempSubRegionSetMap = bcTempStatusMap[regionName][subRegionName];
        if( tempSubRegionSetMap.count( setName ) == 0 )
        {
          bcConsistent = false;
          GEOS_WARNING( BCMessage::missingTemperature( regionName, subRegionName, setName,
                                                       fields::flow::temperature::key() ) );
        }
      }
      if( comp < 0 || comp >= m_numComponents )
      {
        bcConsistent = false;
        GEOS_WARNING( BCMessage::invalidComponentIndex( comp, fs.getName(),
                                                        fields::flow::globalCompFraction::key() ) );
        return; // can't check next part with invalid component id
      }

      ComponentMask< MAX_NC > & compMask = subRegionSetMap[setName];
      if( compMask[comp] )
      {
        bcConsistent = false;
        fsManager.forSubGroups< EquilibriumInitialCondition >( [&] ( EquilibriumInitialCondition const & bc )
        {
          arrayView1d< string const > componentNames = bc.getComponentNames();
          GEOS_WARNING( BCMessage::conflictingComposition( comp, componentNames[comp],
                                                           regionName, subRegionName, setName,
                                                           fields::flow::globalCompFraction::key() ) );
        } );
      }
      compMask.set( comp );
    } );

    // 3.2 Check consistency between composition BC applied to sets
    // Note: for a temperature-only boundary condition, this loop does not do anything
    for( auto const & regionEntry : bcPresCompStatusMap )
    {
      for( auto const & subRegionEntry : regionEntry.second )
      {
        for( auto const & setEntry : subRegionEntry.second )
        {
          ComponentMask< MAX_NC > const & compMask = setEntry.second;

          fsManager.forSubGroups< EquilibriumInitialCondition >( [&] ( EquilibriumInitialCondition const & fs )
          {
            arrayView1d< string const > componentNames = fs.getComponentNames();
            for( int ic = 0; ic < componentNames.size(); ic++ )
            {
              if( !compMask[ic] )
              {
                bcConsistent = false;
                GEOS_WARNING( BCMessage::notAppliedOnRegion( ic, componentNames[ic],
                                                             regionEntry.first, subRegionEntry.first, setEntry.first,
                                                             fields::flow::globalCompFraction::key() ) );
              }
            }
          } );
        }
      }
    }
  } );

  return bcConsistent;
}

void CompositionalMultiphaseBase::applyDirichletBC( real64 const time_n,
                                                    real64 const dt,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  // Only validate BC at the beginning of Newton loop
  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    bool const bcConsistent = validateDirichletBC( domain, time_n + dt );
    GEOS_ERROR_IF( !bcConsistent, GEOS_FMT( "CompositionalMultiphaseBase {}: inconsistent boundary conditions", getDataContext() ) );
  }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {

    // 1. Apply pressure Dirichlet BCs, store in a separate field
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::pressure::key(), fields::flow::bcPressure::key() );
    // 2. Apply composition BC (global component fraction) and store them for constitutive call
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::globalCompFraction::key(), fields::flow::globalCompFraction::key() );
    // 3. Apply temperature Dirichlet BCs, store in a separate field
    if( m_isThermal )
    {
      applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                               fields::flow::temperature::key(), fields::flow::bcTemperature::key() );
    }

    globalIndex const rankOffset = dofManager.rankOffset();
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    // 4. Call constitutive update, back-calculate target global component densities and apply to the system
    fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                             mesh,
                                             fields::flow::pressure::key(),
                                             [&] ( FieldSpecificationBase const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      // in the isothermal case, we use the reservoir temperature to enforce the boundary condition
      // in the thermal case, the validation function guarantees that temperature has been provided
      string const temperatureKey = m_isThermal ? fields::flow::bcTemperature::key() : fields::flow::temperature::key();

      arrayView1d< real64 const > const bcPres =
        subRegion.getReference< array1d< real64 > >( fields::flow::bcPressure::key() );
      arrayView1d< real64 const > const bcTemp =
        subRegion.getReference< array1d< real64 > >( temperatureKey );
      arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::globalCompFraction::key() );

      constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
      {
        using FluidType = TYPEOFREF( castedFluid );
        using ExecPolicy = typename FluidType::exec_policy;
        typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

        thermalCompositionalMultiphaseBaseKernels::
          FluidUpdateKernel::
          launch< ExecPolicy >( targetSet,
                                fluidWrapper,
                                bcPres,
                                bcTemp,
                                compFrac );
      } );

      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< real64 const > const pres =
        subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const compDens =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::globalCompDensity::key() );
      arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();

      integer const numComp = m_numComponents;
      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        // 4.1. Apply pressure value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    bcPres[ei],
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;

        // 4.2. For each component, apply target global density value
        for( integer ic = 0; ic < numComp; ++ic )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      totalDens[ei][0] * compFrac[ei][ic],
                                                      compDens[ei][ic] );
          localRhs[localRow + ic + 1] = rhsValue;
        }
      } );
    } );

    // 5. Apply temperature to the system
    if( m_isThermal )
    {
      fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                               mesh,
                                               fields::flow::temperature::key(),
                                               [&] ( FieldSpecificationBase const &,
                                                     string const &,
                                                     SortedArrayView< localIndex const > const & targetSet,
                                                     ElementSubRegionBase & subRegion,
                                                     string const & )
      {
        arrayView1d< integer const > const ghostRank =
          subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
        arrayView1d< globalIndex const > const dofNumber =
          subRegion.getReference< array1d< globalIndex > >( dofKey );
        arrayView1d< real64 const > const bcTemp =
          subRegion.getReference< array1d< real64 > >( fields::flow::bcTemperature::key() );
        arrayView1d< real64 const > const temp =
          subRegion.getReference< array1d< real64 > >( fields::flow::temperature::key() );

        integer const numComp = m_numComponents;
        forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
        {
          localIndex const ei = targetSet[a];
          if( ghostRank[ei] >= 0 )
          {
            return;
          }

          globalIndex const dofIndex = dofNumber[ei];
          localIndex const localRow = dofIndex - rankOffset;
          real64 rhsValue;

          // 4.2. Apply temperature value to the matrix/rhs
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + numComp + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcTemp[ei],
                                                      temp[ei] );
          localRhs[localRow + numComp + 1] = rhsValue;
        } );
      } );
    }
  } );
}

void CompositionalMultiphaseBase::keepVariablesConstantDuringInitStep( real64 const time,
                                                                       real64 const dt,
                                                                       DofManager const & dofManager,
                                                                       DomainPartition & domain,
                                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                       arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time, dt );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      globalIndex const rankOffset = dofManager.rankOffset();
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::flow::globalCompDensity >();

      integer const numComp = m_numComponents;
      integer const isThermal = m_isThermal;
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        // 4.1. Apply pressure value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    pres[ei], // freeze the current pressure value
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;

        // 4.2. Apply temperature value to the matrix/rhs
        if( isThermal )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + numComp + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      temp[ei], // freeze the current temperature value
                                                      temp[ei] );
          localRhs[localRow + numComp + 1] = rhsValue;
        }

        // 4.3. For each component, apply target global density value
        for( integer ic = 0; ic < numComp; ++ic )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      compDens[ei][ic], // freeze the current component density values
                                                      compDens[ei][ic] );
          localRhs[localRow + ic + 1] = rhsValue;
        }
      } );
    } );
  } );
}

void CompositionalMultiphaseBase::chopNegativeDensities( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  using namespace isothermalCompositionalMultiphaseBaseKernels;

  integer const numComp = m_numComponents;
  real64 const minCompDens = m_minCompDens;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      arrayView2d< real64, compflow::USD_COMP > const compDens =
        subRegion.getField< fields::flow::globalCompDensity >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          for( integer ic = 0; ic < numComp; ++ic )
          {
            if( compDens[ei][ic] < minCompDens )
            {
              compDens[ei][ic] = minCompDens;
            }
          }
        }
      } );
    } );
  } );
}

real64 CompositionalMultiphaseBase::setNextDtBasedOnStateChange( real64 const & currentDt,
                                                                 DomainPartition & domain )
{
  if( m_targetRelativePresChange >= 1.0 &&
      m_targetPhaseVolFracChange >= 1.0 &&
      m_targetRelativeCompDensChange >= 1.0 &&
      ( !m_isThermal || m_targetRelativeTempChange >= 1.0 ) )
  {
    return LvArray::NumericLimits< real64 >::max;
  }

  real64 maxRelativePresChange = 0.0;
  real64 maxRelativeTempChange = 0.0;
  real64 maxAbsolutePhaseVolFracChange = 0.0;
  real64 maxRelativeCompDensChange = 0.0;

  integer const numPhase = m_numPhases;
  integer const numComp = m_numComponents;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const pres_n = subRegion.getField< fields::flow::pressure_n >();
      arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
      arrayView1d< real64 const > const temp_n = subRegion.getField< fields::flow::temperature_n >();
      arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::flow::phaseVolumeFraction >();
      arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac_n =
        subRegion.getField< fields::flow::phaseVolumeFraction_n >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens =
        subRegion.getField< fields::flow::globalCompDensity >();
      arrayView2d< real64, compflow::USD_COMP > const compDens_n =
        subRegion.getField< fields::flow::globalCompDensity_n >();

      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxPresChange( 0.0 );
      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxTempChange( 0.0 );
      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxPhaseVolFracChange( 0.0 );
      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxCompDensChange( 0.0 );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          // switch from relative to absolute when values less than 1
          subRegionMaxPresChange.max( LvArray::math::abs( pres[ei] - pres_n[ei] ) / LvArray::math::max( LvArray::math::abs( pres_n[ei] ), 1.0 ) );
          subRegionMaxTempChange.max( LvArray::math::abs( temp[ei] - temp_n[ei] ) / LvArray::math::max( LvArray::math::abs( temp_n[ei] ), 1.0 ) );
          for( integer ip = 0; ip < numPhase; ++ip )
          {
            subRegionMaxPhaseVolFracChange.max( LvArray::math::abs( phaseVolFrac[ei][ip] - phaseVolFrac_n[ei][ip] ) );
          }
          for( integer ic = 0; ic < numComp; ++ic )
          {
            subRegionMaxCompDensChange.max( LvArray::math::abs( compDens[ei][ic] - compDens_n[ei][ic] ) / LvArray::math::max( LvArray::math::abs( compDens_n[ei][ic] ), 1.0 ) );
          }
        }
      } );

      maxRelativePresChange = LvArray::math::max( maxRelativePresChange, subRegionMaxPresChange.get() );
      maxRelativeTempChange = LvArray::math::max( maxRelativeTempChange, subRegionMaxTempChange.get() );
      maxAbsolutePhaseVolFracChange = LvArray::math::max( maxAbsolutePhaseVolFracChange, subRegionMaxPhaseVolFracChange.get() );
      maxRelativeCompDensChange = LvArray::math::max( maxRelativeCompDensChange, subRegionMaxCompDensChange.get() );

    } );
  } );

  maxRelativePresChange = MpiWrapper::max( maxRelativePresChange );
  maxAbsolutePhaseVolFracChange = MpiWrapper::max( maxAbsolutePhaseVolFracChange );

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: max relative pressure change during time step = {} %",
                                                           getName(), GEOS_FMT( "{:.{}f}", 100*maxRelativePresChange, 3 ) ) );
  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: max absolute phase volume fraction change during time step = {}",
                                                           getName(), GEOS_FMT( "{:.{}f}", maxAbsolutePhaseVolFracChange, 3 ) ) );

  if( m_targetRelativeCompDensChange < LvArray::NumericLimits< real64 >::max )
  {
    maxRelativeCompDensChange = MpiWrapper::max( maxRelativeCompDensChange );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: max relative component density change during time step = {} %",
                                                             getName(), GEOS_FMT( "{:.{}f}", 100*maxRelativeCompDensChange, 3 ) ) );
  }

  if( m_isThermal )
  {
    maxRelativeTempChange = MpiWrapper::max( maxRelativeTempChange );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: max relative temperature change during time step = {} %",
                                                             getName(), GEOS_FMT( "{:.{}f}", 100*maxRelativeTempChange, 3 ) ) );
  }

  real64 const eps = LvArray::NumericLimits< real64 >::epsilon;

  real64 const nextDtPressure = currentDt * ( 1.0 + m_solutionChangeScalingFactor ) * m_targetRelativePresChange
                                / std::max( eps, maxRelativePresChange + m_solutionChangeScalingFactor * m_targetRelativePresChange );
  if( m_nonlinearSolverParameters.getLogLevel() > 0 )
    GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on pressure change = {}", getName(), nextDtPressure ));
  real64 const nextDtPhaseVolFrac = currentDt * ( 1.0 + m_solutionChangeScalingFactor ) * m_targetPhaseVolFracChange
                                    / std::max( eps, maxAbsolutePhaseVolFracChange + m_solutionChangeScalingFactor * m_targetPhaseVolFracChange );
  if( m_nonlinearSolverParameters.getLogLevel() > 0 )
    GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on phase volume fraction change = {}", getName(), nextDtPhaseVolFrac ));
  real64 nextDtCompDens = LvArray::NumericLimits< real64 >::max;
  if( m_targetRelativeCompDensChange < LvArray::NumericLimits< real64 >::max )
  {
    nextDtCompDens = currentDt * ( 1.0 + m_solutionChangeScalingFactor ) * m_targetRelativeCompDensChange
                     / std::max( eps, maxRelativeCompDensChange + m_solutionChangeScalingFactor * m_targetRelativeCompDensChange );
    if( m_nonlinearSolverParameters.getLogLevel() > 0 )
      GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on component density change = {}", getName(), nextDtCompDens ));
  }
  real64 nextDtTemperature = LvArray::NumericLimits< real64 >::max;
  if( m_isThermal )
  {
    nextDtTemperature = currentDt * ( 1.0 + m_solutionChangeScalingFactor ) * m_targetRelativeTempChange
                        / std::max( eps, maxRelativeTempChange + m_solutionChangeScalingFactor * m_targetRelativeTempChange );
    if( m_nonlinearSolverParameters.getLogLevel() > 0 )
      GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on temperature change = {}", getName(), nextDtPhaseVolFrac ));
  }

  return std::min( std::min( nextDtPressure, std::min( nextDtPhaseVolFrac, nextDtCompDens ) ), nextDtTemperature );
}

real64 CompositionalMultiphaseBase::setNextDtBasedOnCFL( const geos::real64 & currentDt, geos::DomainPartition & domain )
{

  real64 maxPhaseCFL, maxCompCFL;

  computeCFLNumbers( domain, currentDt, maxPhaseCFL, maxCompCFL );

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: max phase CFL number = {}", getName(), maxPhaseCFL ) );
  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: max component CFL number = {} ", getName(), maxCompCFL ) );

  return std::min( m_targetFlowCFL * currentDt / maxCompCFL,
                   m_targetFlowCFL * currentDt / maxPhaseCFL );
}

void CompositionalMultiphaseBase::computeCFLNumbers( geos::DomainPartition & domain, const geos::real64 & dt,
                                                     geos::real64 & maxPhaseCFL, geos::real64 & maxCompCFL )
{
  GEOS_MARK_FUNCTION;

  integer const numPhases = numFluidPhases();
  integer const numComps = numFluidComponents();

  // Step 1: reset the arrays involved in the computation of CFL numbers
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView2d< real64, compflow::USD_PHASE > const & phaseOutflux =
        subRegion.getField< fields::flow::phaseOutflux >();
      arrayView2d< real64, compflow::USD_COMP > const & compOutflux =
        subRegion.getField< fields::flow::componentOutflux >();
      phaseOutflux.zero();
      compOutflux.zero();
    } );

    // Step 2: compute the total volumetric outflux for each reservoir cell by looping over faces
    NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( getDiscretizationName() );

    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::PermeabilityAccessors permeabilityAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      CFLFluxKernel::RelPermAccessors relPermAccessors( mesh.getElemManager(), getName() );

    // TODO: find a way to compile with this modifiable accessors in CompFlowAccessors, and remove them from here
    ElementRegionManager::ElementViewAccessor< arrayView2d< real64, compflow::USD_PHASE > > const phaseOutfluxAccessor =
      mesh.getElemManager().constructViewAccessor< array2d< real64, compflow::LAYOUT_PHASE >,
                                                   arrayView2d< real64, compflow::USD_PHASE > >( fields::flow::phaseOutflux::key() );

    ElementRegionManager::ElementViewAccessor< arrayView2d< real64, compflow::USD_COMP > > const compOutfluxAccessor =
      mesh.getElemManager().constructViewAccessor< array2d< real64, compflow::LAYOUT_COMP >,
                                                   arrayView2d< real64, compflow::USD_COMP > >( fields::flow::componentOutflux::key() );


    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      // While this kernel is waiting for a factory class, pass all the accessors here
      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector1
      < isothermalCompositionalMultiphaseFVMKernels::CFLFluxKernel >( numComps,
                                                                      numPhases,
                                                                      dt,
                                                                      stencilWrapper,
                                                                      compFlowAccessors.get( fields::flow::pressure{} ),
                                                                      compFlowAccessors.get( fields::flow::gravityCoefficient{} ),
                                                                      compFlowAccessors.get( fields::flow::phaseVolumeFraction{} ),
                                                                      permeabilityAccessors.get( fields::permeability::permeability{} ),
                                                                      permeabilityAccessors.get( fields::permeability::dPerm_dPressure{} ),
                                                                      relPermAccessors.get( fields::relperm::phaseRelPerm{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseViscosity{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseDensity{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseMassDensity{} ),
                                                                      multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ),
                                                                      phaseOutfluxAccessor.toNestedView(),
                                                                      compOutfluxAccessor.toNestedView() );
    } );
  } );

  // Step 3: finalize the (cell-based) computation of the CFL numbers
  real64 localMaxPhaseCFLNumber = 0.0;
  real64 localMaxCompCFLNumber = 0.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux =
        subRegion.getField< fields::flow::phaseOutflux >();
      arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux =
        subRegion.getField< fields::flow::componentOutflux >();

      arrayView1d< real64 > const & phaseCFLNumber = subRegion.getField< fields::flow::phaseCFLNumber >();
      arrayView1d< real64 > const & compCFLNumber = subRegion.getField< fields::flow::componentCFLNumber >();

      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

      arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
        subRegion.getField< fields::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
        subRegion.getField< fields::flow::globalCompFraction >();
      arrayView2d< real64, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::flow::phaseVolumeFraction >();

      Group const & constitutiveModels = subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

      string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = constitutiveModels.getGroup< MultiFluidBase >( fluidName );
      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc = fluid.phaseViscosity();

      string const & relpermName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relperm = constitutiveModels.getGroup< RelativePermeabilityBase >( relpermName );
      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm = relperm.phaseRelPerm();
      arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac = relperm.dPhaseRelPerm_dPhaseVolFraction();

      string const & solidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = constitutiveModels.getGroup< CoupledSolidBase >( solidName );
      arrayView2d< real64 const > const & porosity = solid.getPorosity();

      real64 subRegionMaxPhaseCFLNumber = 0.0;
      real64 subRegionMaxCompCFLNumber = 0.0;

      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector2
      < isothermalCompositionalMultiphaseFVMKernels::CFLKernel >( numComps, numPhases,
                                                                  subRegion.size(),
                                                                  volume,
                                                                  porosity,
                                                                  compDens,
                                                                  compFrac,
                                                                  phaseVolFrac,
                                                                  phaseRelPerm,
                                                                  dPhaseRelPerm_dPhaseVolFrac,
                                                                  phaseVisc,
                                                                  phaseOutflux,
                                                                  compOutflux,
                                                                  phaseCFLNumber,
                                                                  compCFLNumber,
                                                                  subRegionMaxPhaseCFLNumber,
                                                                  subRegionMaxCompCFLNumber );

      localMaxPhaseCFLNumber = LvArray::math::max( localMaxPhaseCFLNumber, subRegionMaxPhaseCFLNumber );
      localMaxCompCFLNumber = LvArray::math::max( localMaxCompCFLNumber, subRegionMaxCompCFLNumber );

    } );
  } );

  maxPhaseCFL = MpiWrapper::max( localMaxPhaseCFLNumber );
  maxCompCFL = MpiWrapper::max( localMaxCompCFLNumber );

}


void CompositionalMultiphaseBase::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      arrayView1d< real64 > const & pres =
        subRegion.template getField< fields::flow::pressure >();
      arrayView1d< real64 const > const & pres_n =
        subRegion.template getField< fields::flow::pressure_n >();
      pres.setValues< parallelDevicePolicy<> >( pres_n );

      arrayView2d< real64, compflow::USD_COMP > const & compDens =
        subRegion.template getField< fields::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const & compDens_n =
        subRegion.template getField< fields::flow::globalCompDensity_n >();
      compDens.setValues< parallelDevicePolicy<> >( compDens_n );

      if( m_isThermal )
      {
        arrayView1d< real64 > const & temp =
          subRegion.template getField< fields::flow::temperature >();
        arrayView1d< real64 const > const & temp_n =
          subRegion.template getField< fields::flow::temperature_n >();
        temp.setValues< parallelDevicePolicy<> >( temp_n );
      }

      // update porosity, permeability
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updateFluidState( subRegion );
      // for thermal simulations, update solid internal energy
      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
      }

    } );
  } );
}

void CompositionalMultiphaseBase::implicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  // Step 1: save the converged aquifer state
  // note: we have to save the aquifer state **before** updating the pressure,
  // otherwise the aquifer flux is saved with the wrong pressure time level
  saveAquiferConvergedState( time, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // update deltaPressure
      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const initPres = subRegion.getField< fields::flow::initialPressure >();
      arrayView1d< real64 > const deltaPres = subRegion.getField< fields::flow::deltaPressure >();
      isothermalCompositionalMultiphaseBaseKernels::StatisticsKernel::
        saveDeltaPressure< parallelDevicePolicy<> >( subRegion.size(), pres, initPres, deltaPres );

      // Step 2: save the converged fluid state
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluidMaterial = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      fluidMaterial.saveConvergedState();

      // Step 3: save the converged solid state
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
      if( m_keepVariablesConstantDuringInitStep )
      {
        porousMaterial.ignoreConvergedState(); // newPorosity <- porosity_n
      }
      else
      {
        porousMaterial.saveConvergedState(); // porosity_n <- porosity
      }

      // Step 4: save converged state for the relperm model to handle hysteresis
      arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::flow::phaseVolumeFraction >();
      string const & relPermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relPermMaterial =
        getConstitutiveModel< RelativePermeabilityBase >( subRegion, relPermName );
      relPermMaterial.saveConvergedPhaseVolFractionState( phaseVolFrac );

      // Step 5: if capillary pressure is supported, send the converged porosity and permeability to the capillary pressure model
      // note: this is needed when the capillary pressure depends on porosity and permeability (Leverett J-function for instance)
      if( m_hasCapPressure )
      {
        arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

        string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
        PermeabilityBase const & permeabilityMaterial =
          getConstitutiveModel< PermeabilityBase >( subRegion, permName );
        arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

        string const & capPressName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        CapillaryPressureBase const & capPressureMaterial =
          getConstitutiveModel< CapillaryPressureBase >( subRegion, capPressName );
        capPressureMaterial.saveConvergedRockState( porosity, permeability );
      }

      // Step 6: if the thermal option is on, send the converged porosity and phase volume fraction to the thermal conductivity model
      // note: this is needed because the phaseVolFrac-weighted thermal conductivity treats phaseVolumeFraction explicitly for now
      if( m_isThermal )
      {
        arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

        string const & thermName = subRegion.getReference< string >( viewKeyStruct::thermalConductivityNamesString() );
        MultiPhaseThermalConductivityBase const & thermalConductivityMaterial =
          getConstitutiveModel< MultiPhaseThermalConductivityBase >( subRegion, thermName );
        thermalConductivityMaterial.saveConvergedRockFluidState( porosity, phaseVolFrac );
      }

      // Step 7: if the diffusion and/or dispersion is/are supported, update the two models explicity
      if( m_hasDiffusion )
      {
        string const & diffusionName = subRegion.getReference< string >( viewKeyStruct::diffusionNamesString() );
        DiffusionBase const & diffusionMaterial = getConstitutiveModel< DiffusionBase >( subRegion, diffusionName );
        arrayView1d< real64 const > const temperature = subRegion.template getField< fields::flow::temperature >();
        diffusionMaterial.saveConvergedTemperatureState( temperature );
      }
      if( m_hasDispersion )
      {
        string const & dispersionName = subRegion.getReference< string >( viewKeyStruct::dispersionNamesString() );
        DispersionBase const & dispersionMaterial = getConstitutiveModel< DispersionBase >( subRegion, dispersionName );
        GEOS_UNUSED_VAR( dispersionMaterial );
        // TODO: compute the total velocity here
        //dispersionMaterial.saveConvergedVelocitySate( totalVelovity );
      }
    } );
  } );
}

void CompositionalMultiphaseBase::saveConvergedState( ElementSubRegionBase & subRegion ) const
{
  FlowSolverBase::saveConvergedState( subRegion );

  arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
    subRegion.template getField< fields::flow::globalCompDensity >();
  arrayView2d< real64, compflow::USD_COMP > const & compDens_n =
    subRegion.template getField< fields::flow::globalCompDensity_n >();
  compDens_n.setValues< parallelDevicePolicy<> >( compDens );

  arrayView2d< real64 const, compflow::USD_COMP > const & compAmount =
    subRegion.template getField< fields::flow::compAmount >();
  arrayView2d< real64, compflow::USD_COMP > const & compAmount_n =
    subRegion.template getField< fields::flow::compAmount_n >();
  compAmount_n.setValues< parallelDevicePolicy<> >( compAmount );

  if( m_isFixedStressPoromechanicsUpdate )
  {
    arrayView2d< real64, compflow::USD_COMP > const & compDens_k =
      subRegion.template getField< fields::flow::globalCompDensity_k >();
    compDens_k.setValues< parallelDevicePolicy<> >( compDens );
  }
}

void CompositionalMultiphaseBase::saveSequentialIterationState( DomainPartition & domain )
{
  FlowSolverBase::saveSequentialIterationState( domain );

  integer const numComp = m_numComponents;

  real64 maxCompDensChange = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      arrayView2d< real64 const, compflow::USD_COMP >
      const compDens = subRegion.getField< fields::flow::globalCompDensity >();
      arrayView2d< real64, compflow::USD_COMP >
      const compDens_k = subRegion.getField< fields::flow::globalCompDensity_k >();

      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxCompDensChange( 0.0 );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=]
                                        GEOS_HOST_DEVICE ( localIndex
                                                           const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          for( integer ic = 0; ic < numComp; ++ic )
          {
            subRegionMaxCompDensChange.max( LvArray::math::abs( compDens[ei][ic] - compDens_k[ei][ic] ) );
            compDens_k[ei][ic] = compDens[ei][ic];
          }
        }
      } );

      maxCompDensChange = LvArray::math::max( maxCompDensChange, subRegionMaxCompDensChange.get() );
    } );
  } );

  m_sequentialCompDensChange = MpiWrapper::max( maxCompDensChange ); // store to be later used for convergence check
}

void CompositionalMultiphaseBase::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  real64 maxDeltaPhaseVolFrac = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                             auto & subRegion )
    {
      // update porosity, permeability, and solid internal energy
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      real64 const deltaPhaseVolFrac = updateFluidState( subRegion );
      maxDeltaPhaseVolFrac = LvArray::math::max( maxDeltaPhaseVolFrac, deltaPhaseVolFrac );
      // for thermal, update solid internal energy
      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
        updateEnergy( subRegion );
      }
    } );
  } );

  maxDeltaPhaseVolFrac = MpiWrapper::max( maxDeltaPhaseVolFrac );

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution,
                              GEOS_FMT( "        {}: Max phase volume fraction change = {}", getName(), fmt::format( "{:.{}f}", maxDeltaPhaseVolFrac, 4 ) ) );
}

bool CompositionalMultiphaseBase::checkSequentialSolutionIncrements( DomainPartition & domain ) const
{
  bool isConverged = FlowSolverBase::checkSequentialSolutionIncrements( domain );

  string const unit = m_useMass ? "kg/m3" : "mol/m3";
  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Convergence, GEOS_FMT( "    {}: Max component density change during outer iteration: {} {}",
                                                              getName(), fmt::format( "{:.{}f}", m_sequentialCompDensChange, 3 ), unit ) );

  return isConverged && (m_sequentialCompDensChange < m_maxSequentialCompDensChange);
}

real64 CompositionalMultiphaseBase::setNextDt( const geos::real64 & currentDt, geos::DomainPartition & domain )
{

  if( m_targetFlowCFL < 0 )
    return PhysicsSolverBase::setNextDt( currentDt, domain );
  else
    return setNextDtBasedOnCFL( currentDt, domain );
}

} // namespace geos
