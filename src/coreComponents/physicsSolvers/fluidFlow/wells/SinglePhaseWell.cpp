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
 * @file SinglePhaseWell.cpp
 */

#include "SinglePhaseWell.hpp"

#include "common/DataTypes.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidSelector.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluidUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/SolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace singlePhaseWellKernels;

SinglePhaseWell::SinglePhaseWell( const string & name,
                                  Group * const parent ):
  WellSolverBase( name, parent )
{
  m_numDofPerWellElement = 2;
  m_numDofPerResElement = 1;
  m_numPhases = 1;
  m_numComponents = 1;

  this->registerWrapper( FlowSolverBase::viewKeyStruct::allowNegativePressureString(), &m_allowNegativePressure ).
    setApplyDefaultValue( 1 ). // negative pressure is allowed by default
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating if negative pressure is allowed" );
}

void SinglePhaseWell::registerDataOnMesh( Group & meshBodies )
{
  WellSolverBase::registerDataOnMesh( meshBodies );

  // loop over the wells
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      fluidName = getConstitutiveName< SingleFluidBase >( subRegion );
      GEOS_ERROR_IF( fluidName.empty(), GEOS_FMT( "{}: Fluid model not found on subregion {}",
                                                  getDataContext(), subRegion.getName() ) );

      subRegion.registerField< fields::well::connectionRate_n >( getName() );
      subRegion.registerField< fields::well::connectionRate >( getName() );

      PerforationData & perforationData = *subRegion.getPerforationData();
      perforationData.registerField< fields::well::perforationRate >( getName() );
      perforationData.registerField< fields::well::dPerforationRate_dPres >( getName() ).
        reference().resizeDimension< 1 >( 2 );

      WellControls & wellControls = getWellControls( subRegion );
      wellControls.registerWrapper< real64 >( viewKeyStruct::currentBHPString() );
      wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentBHP_dPresString() );

      wellControls.registerWrapper< real64 >( viewKeyStruct::currentVolRateString() );
      wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentVolRate_dPresString() );
      wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentVolRate_dRateString() );

      // write rates output header
      if( m_writeCSV > 0 && subRegion.isLocallyOwned())
      {
        string const wellControlsName = wellControls.getName();
        integer const useSurfaceConditions = wellControls.useSurfaceConditions();
        string const conditionKey = useSurfaceConditions ? "surface" : "reservoir";
        string const unitKey = useSurfaceConditions ? "s" : "r";
        // format: time,bhp,total_rate,total_vol_rate
        std::ofstream outputFile( m_ratesOutputDir + "/" + wellControlsName + ".csv" );
        outputFile << "Time [s],BHP [Pa],Total rate [kg/s],Total " << conditionKey << " volumetric rate ["<<unitKey<<"m3/s]" << std::endl;
        outputFile.close();
      }
    } );
  } );
}

string SinglePhaseWell::resElementDofName() const
{
  return SinglePhaseBase::viewKeyStruct::elemDofFieldString();
}

void SinglePhaseWell::validateWellConstraints( real64 const & time_n,
                                               real64 const & dt,
                                               WellElementSubRegion const & subRegion )
{
  WellControls const & wellControls = getWellControls( subRegion );
  WellControls::Control const currentControl = wellControls.getControl();
  real64 const targetTotalRate = wellControls.getTargetTotalRate( time_n + dt );
  real64 const targetPhaseRate = wellControls.getTargetPhaseRate( time_n + dt );
  GEOS_THROW_IF( currentControl == WellControls::Control::PHASEVOLRATE,
                 "WellControls " << wellControls.getDataContext() <<
                 ": Phase rate control is not available for SinglePhaseWell",
                 InputError );
  // The user always provides positive rates, but these rates are later multiplied by -1 internally for producers
  GEOS_THROW_IF( ( ( wellControls.isInjector() && targetTotalRate < 0.0 ) ||
                   ( wellControls.isProducer() && targetTotalRate > 0.0) ),
                 "WellControls " << wellControls.getDataContext() <<
                 ": Target total rate cannot be negative",
                 InputError );
  GEOS_THROW_IF( !isZero( targetPhaseRate ),
                 "WellControls " << wellControls.getDataContext() <<
                 ": Target phase rate cannot be used for SinglePhaseWell",
                 InputError );
}

void SinglePhaseWell::updateBHPForConstraint( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const pres =
    subRegion.getField< fields::well::pressure >();

  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getField< fields::well::gravityCoefficient >();

  // fluid data
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
  arrayView2d< real64 const > const & dens = fluid.density();
  arrayView2d< real64 const > const & dDens_dPres = fluid.dDensity_dPressure();

  // control data

  WellControls & wellControls = getWellControls( subRegion );
  string const wellControlsName = wellControls.getName();
  integer const logLevel = wellControls.getLogLevel();
  real64 const & refGravCoef = wellControls.getReferenceGravityCoef();

  real64 & currentBHP =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );
  real64 & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentBHP_dPresString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [pres,
                              dens,
                              dDens_dPres,
                              wellElemGravCoef,
                              &currentBHP,
                              &dCurrentBHP_dPres,
                              &iwelemRef,
                              &refGravCoef] ( localIndex const )
  {
    currentBHP = pres[iwelemRef] + dens[iwelemRef][0] * ( refGravCoef - wellElemGravCoef[iwelemRef] );
    dCurrentBHP_dPres = 1.0 + dDens_dPres[iwelemRef][0] * ( refGravCoef - wellElemGravCoef[iwelemRef] );
  } );

  if( logLevel >= 2 )
  {
    GEOS_LOG_RANK( GEOS_FMT( "{}: The BHP (at the specified reference elevation) = {} Pa",
                             wellControlsName, currentBHP ) );
  }

}

void SinglePhaseWell::updateVolRateForConstraint( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const pres =
    subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const & connRate =
    subRegion.getField< fields::well::connectionRate >();

  // fluid data

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
  arrayView2d< real64 const > const & dens = fluid.density();
  arrayView2d< real64 const > const & dDens_dPres = fluid.dDensity_dPressure();

  // control data

  WellControls & wellControls = getWellControls( subRegion );
  string const wellControlsName = wellControls.getName();
  integer const logLevel = wellControls.getLogLevel();
  integer const useSurfaceConditions = wellControls.useSurfaceConditions();
  real64 const & surfacePres = wellControls.getSurfacePressure();

  real64 & currentVolRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString() );
  real64 & dCurrentVolRate_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dPresString() );
  real64 & dCurrentVolRate_dRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dRateString() );

  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [fluidWrapper,
                                pres,
                                connRate,
                                dens,
                                dDens_dPres,
                                &useSurfaceConditions,
                                &surfacePres,
                                &currentVolRate,
                                &dCurrentVolRate_dPres,
                                &dCurrentVolRate_dRate,
                                &iwelemRef,
                                &logLevel,
                                &wellControlsName] ( localIndex const )
    {
      //    We need to evaluate the density as follows:
      //      - Surface conditions: using the surface pressure provided by the user
      //      - Reservoir conditions: using the pressure in the top element

      if( useSurfaceConditions )
      {
        // we need to compute the surface density
        fluidWrapper.update( iwelemRef, 0, surfacePres );
        if( logLevel >= 2 )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: surface density computed with P_surface = {} Pa",
                                   wellControlsName, surfacePres ) );
#ifdef GEOS_USE_HIP
          GEOS_UNUSED_VAR( wellControlsName );
#endif
        }
      }
      else
      {
        real64 const refPres = pres[iwelemRef];
        fluidWrapper.update( iwelemRef, 0, refPres );
      }

      real64 const densInv = 1.0 / dens[iwelemRef][0];
      currentVolRate = connRate[iwelemRef] * densInv;
      dCurrentVolRate_dPres = -( useSurfaceConditions ==  0 ) * dDens_dPres[iwelemRef][0] * currentVolRate * densInv;
      dCurrentVolRate_dRate = densInv;

      if( logLevel >= 2 && useSurfaceConditions )
      {
        GEOS_LOG_RANK( GEOS_FMT( "{}: total fluid density at surface conditions = {} kg/sm3, total rate = {} kg/s, total surface volumetric rate = {} sm3/s",
                                 wellControlsName, dens[iwelemRef][0], connRate[iwelemRef], currentVolRate ) );
      }
    } );
  } );
}

void SinglePhaseWell::updateFluidModel( WellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const temp = subRegion.getField< fields::well::temperature >();

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );

  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    singlePhaseBaseKernels::FluidUpdateKernel::launch( fluidWrapper, pres, temp );
  } );
}

real64 SinglePhaseWell::updateSubRegionState( WellElementSubRegion & subRegion )
{
  // update volumetric rates for the well constraints
  // Warning! This must be called before updating the fluid model
  updateVolRateForConstraint( subRegion );

  // update density in the well elements
  updateFluidModel( subRegion );

  // update the current BHP
  updateBHPForConstraint( subRegion );

  // note: the perforation rates are updated separately
  return 0.0;  // change in phasevolume fraction doesnt apply
}

void SinglePhaseWell::initializeWells( DomainPartition & domain, real64 const & time_n, real64 const & dt )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );
  // different functionality than for compositional
  // compositional has better treatment of well initialization in cases where well starts later in sim
  // logic will be incorporated into single phase in different push
  if( time_n > 0 )
  {
    return;
  }
  // loop over the wells
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      WellControls const & wellControls = getWellControls( subRegion );
      PerforationData const & perforationData = *subRegion.getPerforationData();

      // get the info stored on well elements
      arrayView1d< real64 const > const wellElemGravCoef =
        subRegion.getField< fields::well::gravityCoefficient >();

      // get well primary variables on well elements
      arrayView1d< real64 > const wellElemPressure =
        subRegion.getField< fields::well::pressure >();
      arrayView1d< real64 > const connRate =
        subRegion.getField< fields::well::connectionRate >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const resElementRegion =
        perforationData.getField< fields::perforation::reservoirElementRegion >();
      arrayView1d< localIndex const > const resElementSubRegion =
        perforationData.getField< fields::perforation::reservoirElementSubRegion >();
      arrayView1d< localIndex const > const resElementIndex =
        perforationData.getField< fields::perforation::reservoirElementIndex >();

      arrayView1d< real64 const > const & perfGravCoef =
        perforationData.getField< fields::well::gravityCoefficient >();

      // TODO: change the way we access the flowSolver here
      SinglePhaseBase const & flowSolver = getParent().getGroup< SinglePhaseBase >( getFlowSolverName() );
      PresInitializationKernel::SinglePhaseFlowAccessors resSinglePhaseFlowAccessors( meshLevel.getElemManager(), flowSolver.getName() );
      PresInitializationKernel::SingleFluidAccessors resSingleFluidAccessors( meshLevel.getElemManager(), flowSolver.getName() );

      // 1) Loop over all perforations to compute an average density
      // 2) Initialize the reference pressure
      // 3) Estimate the pressures in the well elements using the average density
      PresInitializationKernel::
        launch( perforationData.size(),
                subRegion.size(),
                perforationData.getNumPerforationsGlobal(),
                wellControls,
                0.0,   // initialization done at t = 0
                resSinglePhaseFlowAccessors.get( fields::flow::pressure{} ),
                resSingleFluidAccessors.get( fields::singlefluid::density{} ),
                resElementRegion,
                resElementSubRegion,
                resElementIndex,
                perfGravCoef,
                wellElemGravCoef,
                wellElemPressure );

      // 4) Recompute the pressure-dependent properties
      // Note: I am leaving that here because I would like to use the perforationRates (computed in UpdateState)
      //       to better initialize the rates
      updateSubRegionState( subRegion );

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
      arrayView2d< real64 const > const & wellElemDens = fluid.density();

      // 5) Estimate the well rates
      RateInitializationKernel::launch( subRegion.size(),
                                        wellControls,
                                        0.0,   // initialization done at t = 0
                                        wellElemDens,
                                        connRate );

    } );

  } );
}

void SinglePhaseWell::assembleFluxTerms( real64 const & time_n,
                                         real64 const & dt,
                                         DomainPartition & domain,
                                         DofManager const & dofManager,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );

  // loop over the wells
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {
      // get a reference to the degree-of-freedom numbers
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< localIndex const > const & nextWellElemIndex =
        subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const connRate =
        subRegion.getField< fields::well::connectionRate >();

      FluxKernel::launch( subRegion.size(),
                          dofManager.rankOffset(),
                          wellElemDofNumber,
                          nextWellElemIndex,
                          connRate,
                          dt,
                          localMatrix,
                          localRhs );
    } );

  } );
}

void SinglePhaseWell::assemblePressureRelations( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition const & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();


    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      WellControls & wellControls = getWellControls( subRegion );

      // get the degrees of freedom numbers, depth, next well elem index
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< real64 const > const & wellElemGravCoef =
        subRegion.getField< fields::well::gravityCoefficient >();
      arrayView1d< localIndex const > const & nextWellElemIndex =
        subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

      // get primary variables on well elements
      arrayView1d< real64 const > const & wellElemPressure =
        subRegion.getField< fields::well::pressure >();

      // get well constitutive data
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
      arrayView2d< real64 const > const & wellElemDensity = fluid.density();
      arrayView2d< real64 const > const & dWellElemDensity_dPres = fluid.dDensity_dPressure();

      localIndex const controlHasSwitched =
        PressureRelationKernel::launch( subRegion.size(),
                                        dofManager.rankOffset(),
                                        subRegion.isLocallyOwned(),
                                        subRegion.getTopWellElementIndex(),
                                        wellControls,
                                        time_n + dt,   // controls evaluated with BHP/rate of the end of the time interval
                                        wellElemDofNumber,
                                        wellElemGravCoef,
                                        nextWellElemIndex,
                                        wellElemPressure,
                                        wellElemDensity,
                                        dWellElemDensity_dPres,
                                        localMatrix,
                                        localRhs );

      if( controlHasSwitched == 1 )
      {
        // Note: if BHP control is not viable, we switch to TOTALVOLRATE
        //       if TOTALVOLRATE is not viable, we switch to BHP

        real64 const timeAtEndOfStep = time_n + dt;

        if( wellControls.getControl() == WellControls::Control::BHP )
        {
          wellControls.switchToTotalRateControl( wellControls.getTargetTotalRate( timeAtEndOfStep ) );
          GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::WellControl, GEOS_FMT( "Control switch for well {} from BHP constraint to rate constraint", subRegion.getName()) );
        }
        else
        {
          wellControls.switchToBHPControl( wellControls.getTargetBHP( timeAtEndOfStep ) );
          GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::WellControl, GEOS_FMT( "Control switch for well {} from rate constraint to BHP constraint", subRegion.getName()) );
        }
      }

    } );
  } );
}

void SinglePhaseWell::assembleAccumulationTerms( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      // get a reference to the degree-of-freedom numbers
      string const wellElemDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellElemDofKey );
      arrayView1d< integer const > const wellElemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const wellElemVolume = subRegion.getElementVolume();

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );
      arrayView2d< real64 const > const wellElemDensity = fluid.density();
      arrayView2d< real64 const > const dWellElemDensity_dPres = fluid.dDensity_dPressure();
      arrayView2d< real64 const > const wellElemDensity_n = fluid.density_n();

      AccumulationKernel::launch( subRegion.size(),
                                  dofManager.rankOffset(),
                                  wellElemDofNumber,
                                  wellElemGhostRank,
                                  wellElemVolume,
                                  wellElemDensity,
                                  dWellElemDensity_dPres,
                                  wellElemDensity_n,
                                  localMatrix,
                                  localRhs );

    } );
  } );
  // then assemble the volume balance equations
  assembleVolumeBalanceTerms( domain, dofManager, localMatrix, localRhs );
}

void SinglePhaseWell::assembleVolumeBalanceTerms( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                  DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                  CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                  arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  // not implemented for single phase flow
}

void SinglePhaseWell::shutDownWell( real64 const time_n,
                                    real64 const dt,
                                    DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      // if the well is open, we don't have to do anything, so we just return
      WellControls const & wellControls = getWellControls( subRegion );
      if( wellControls.isWellOpen( time_n + dt ) )
      {
        return;
      }

      globalIndex const rankOffset = dofManager.rankOffset();

      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      arrayView1d< real64 const > const pres =
        subRegion.getField< fields::well::pressure >();
      arrayView1d< real64 const > const connRate =
        subRegion.getField< fields::well::connectionRate >();

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

        // 4.2. Apply rate value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex + 1,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    connRate[ei], // freeze the current pressure value
                                                    connRate[ei] );
        localRhs[localRow + 1] = rhsValue;

      } );
    } );
  } );
}


void SinglePhaseWell::computePerforationRates( real64 const & time_n,
                                               real64 const & dt, DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    // TODO: change the way we access the flowSolver here
    SinglePhaseBase const & flowSolver = getParent().getGroup< SinglePhaseBase >( getFlowSolverName() );
    PerforationKernel::SinglePhaseFlowAccessors resSinglePhaseFlowAccessors( mesh.getElemManager(), flowSolver.getName() );
    PerforationKernel::SingleFluidAccessors resSingleFluidAccessors( mesh.getElemManager(), flowSolver.getName() );

    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {

      // get the well data
      PerforationData * const perforationData = subRegion.getPerforationData();

      // get the degrees of freedom and depth
      arrayView1d< real64 const > const wellElemGravCoef =
        subRegion.getField< fields::well::gravityCoefficient >();

      // get well primary variables on well elements
      arrayView1d< real64 const > const wellElemPressure =
        subRegion.getField< fields::well::pressure >();

      // get well constitutive data
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );  arrayView2d< real64 const > const wellElemDensity = fluid.density();
      arrayView2d< real64 const > const dWellElemDensity_dPres = fluid.dDensity_dPressure();
      arrayView2d< real64 const > const wellElemViscosity = fluid.viscosity();
      arrayView2d< real64 const > const dWellElemViscosity_dPres = fluid.dViscosity_dPressure();

      // get well variables on perforations
      arrayView1d< real64 const > const perfGravCoef =
        perforationData->getField< fields::well::gravityCoefficient >();
      arrayView1d< localIndex const > const perfWellElemIndex =
        perforationData->getField< fields::perforation::wellElementIndex >();
      arrayView1d< real64 const > const perfTransmissibility =
        perforationData->getField< fields::perforation::wellTransmissibility >();

      arrayView1d< real64 > const perfRate =
        perforationData->getField< fields::well::perforationRate >();
      arrayView2d< real64 > const dPerfRate_dPres =
        perforationData->getField< fields::well::dPerforationRate_dPres >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const resElementRegion =
        perforationData->getField< fields::perforation::reservoirElementRegion >();
      arrayView1d< localIndex const > const resElementSubRegion =
        perforationData->getField< fields::perforation::reservoirElementSubRegion >();
      arrayView1d< localIndex const > const resElementIndex =
        perforationData->getField< fields::perforation::reservoirElementIndex >();

      PerforationKernel::launch( perforationData->size(),
                                 resSinglePhaseFlowAccessors.get( fields::flow::pressure{} ),
                                 resSingleFluidAccessors.get( fields::singlefluid::density{} ),
                                 resSingleFluidAccessors.get( fields::singlefluid::dDensity_dPressure{} ),
                                 resSingleFluidAccessors.get( fields::singlefluid::viscosity{} ),
                                 resSingleFluidAccessors.get( fields::singlefluid::dViscosity_dPressure{} ),
                                 wellElemGravCoef,
                                 wellElemPressure,
                                 wellElemDensity,
                                 dWellElemDensity_dPres,
                                 wellElemViscosity,
                                 dWellElemViscosity_dPres,
                                 perfGravCoef,
                                 perfWellElemIndex,
                                 perfTransmissibility,
                                 resElementRegion,
                                 resElementSubRegion,
                                 resElementIndex,
                                 perfRate,
                                 dPerfRate_dPres );
    } );
  } );
}


real64
SinglePhaseWell::calculateResidualNorm( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  real64 localResidualNorm = 0.0;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {
      real64 subRegionResidualNorm[1]{};

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );

      WellControls const & wellControls = getWellControls( subRegion );

      // step 1: compute the norm in the subRegion

      ResidualNormKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( rankOffset,
                                                   wellDofKey,
                                                   localRhs,
                                                   subRegion,
                                                   fluid,
                                                   wellControls,
                                                   time_n + dt,
                                                   dt,
                                                   m_nonlinearSolverParameters.m_minNormalizer,
                                                   subRegionResidualNorm );

      // step 2: reduction across meshBodies/regions/subRegions

      if( subRegionResidualNorm[0] > localResidualNorm )
      {
        localResidualNorm = subRegionResidualNorm[0];
      }

    } );
  } );


  // step 3: second reduction across MPI ranks

  real64 const residualNorm = MpiWrapper::max( localResidualNorm );

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    std::cout << GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), residualNorm );
  }
  return residualNorm;
}

bool SinglePhaseWell::checkSystemSolution( DomainPartition & domain,
                                           DofManager const & dofManager,
                                           arrayView1d< real64 const > const & localSolution,
                                           real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );
  integer numNegativePressures = 0;
  real64 minPressure = 0.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )

    {
      globalIndex const rankOffset = dofManager.rankOffset();
      // get the degree of freedom numbers on well elements
      arrayView1d< globalIndex const > const & dofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const & pres =
        subRegion.getField< fields::well::pressure >();

      auto const statistics =
        singlePhaseBaseKernels::SolutionCheckKernel::
          launch< parallelDevicePolicy<> >( localSolution, rankOffset, dofNumber, ghostRank, pres, scalingFactor );

      numNegativePressures += statistics.first;
      minPressure = std::min( minPressure, statistics.second );
    } );
  } );

  numNegativePressures = MpiWrapper::sum( numNegativePressures );

  if( numNegativePressures > 0 )
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Number of negative pressure values: {}, minimum value: {} Pa",
                                        getName(), numNegativePressures, fmt::format( "{:.{}f}", minPressure, 3 ) ) );

  return (m_allowNegativePressure || numNegativePressures == 0) ?  1 : 0;
}

void
SinglePhaseWell::applySystemSolution( DofManager const & dofManager,
                                      arrayView1d< real64 const > const & localSolution,
                                      real64 const scalingFactor,
                                      real64 const dt,
                                      DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );
  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               fields::well::pressure::key(),
                               scalingFactor,
                               { m_numDofPerWellElement, 0, 1 } );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               fields::well::connectionRate::key(),
                               scalingFactor,
                               { m_numDofPerWellElement, 1, m_numDofPerWellElement } );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { fields::well::pressure::key(),
                                       fields::well::connectionRate::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );

}

void SinglePhaseWell::resetStateToBeginningOfStep( DomainPartition & domain )
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      // get a reference to the primary variables on well elements
      arrayView1d< real64 > const & wellElemPressure =
        subRegion.getField< fields::well::pressure >();
      arrayView1d< real64 const > const & wellElemPressure_n =
        subRegion.getField< fields::well::pressure_n >();
      wellElemPressure.setValues< parallelDevicePolicy<> >( wellElemPressure_n );

      arrayView1d< real64 > const & connRate =
        subRegion.getField< fields::well::connectionRate >();
      arrayView1d< real64 const > const & connRate_n =
        subRegion.getField< fields::well::connectionRate_n >();
      connRate.setValues< parallelDevicePolicy<> >( connRate_n );

      updateSubRegionState( subRegion );
    } );
  } );
}


void SinglePhaseWell::implicitStepSetup( real64 const & time,
                                         real64 const & dt,
                                         DomainPartition & domain )
{
  WellSolverBase::implicitStepSetup( time, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      arrayView1d< real64 const > const wellElemPressure = subRegion.getField< fields::well::pressure >();
      arrayView1d< real64 > const wellElemPressure_n = subRegion.getField< fields::well::pressure_n >();
      wellElemPressure_n.setValues< parallelDevicePolicy<> >( wellElemPressure );

      arrayView1d< real64 const > const connRate = subRegion.getField< fields::well::connectionRate >();
      arrayView1d< real64 > const connRate_n = subRegion.getField< fields::well::connectionRate_n >();
      connRate_n.setValues< parallelDevicePolicy<> >( connRate );

      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
      fluid.saveConvergedState();

      validateWellConstraints( time, dt, subRegion );

      updateSubRegionState( subRegion );
    } );
  } );
}

void SinglePhaseWell::implicitStepComplete( real64 const & time_n,
                                            real64 const & dt,
                                            DomainPartition & domain )
{
  WellSolverBase::implicitStepComplete( time_n, dt, domain );

  if( getLogLevel() > 0 )
  {
    printRates( time_n, dt, domain );
  }
}

void SinglePhaseWell::printRates( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {

      // the rank that owns the reference well element is responsible for the calculations below.
      if( !subRegion.isLocallyOwned() )
      {
        return;
      }

      localIndex const iwelemRef = subRegion.getTopWellElementIndex();

      // subRegion data

      arrayView1d< real64 const > const & connRate =
        subRegion.getField< fields::well::connectionRate >();

      // control data

      WellControls const & wellControls = getWellControls( subRegion );
      string const wellControlsName = wellControls.getName();

      // format: time,total_rate,total_vol_rate
      std::ofstream outputFile;
      if( m_writeCSV > 0 )
      {
        outputFile.open( m_ratesOutputDir + "/" + wellControlsName + ".csv", std::ios_base::app );
        outputFile << time_n + dt;
      }

      if( !wellControls.isWellOpen( time_n + dt ) )
      {
        GEOS_LOG( GEOS_FMT( "{}: well is shut", wellControlsName ) );
        if( outputFile.is_open())
        {
          // print all zeros in the rates file
          outputFile << ",0.0,0.0,0.0" << std::endl;
          outputFile.close();
        }
        return;
      }

      integer const useSurfaceConditions = wellControls.useSurfaceConditions();

      real64 const & currentBHP =
        wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );
      real64 const & currentTotalVolRate =
        wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString() );

      // bring everything back to host, capture the scalars by reference
      forAll< serialPolicy >( 1, [&useSurfaceConditions,
                                  &currentBHP,
                                  connRate,
                                  &currentTotalVolRate,
                                  &iwelemRef,
                                  &wellControlsName,
                                  &outputFile] ( localIndex const )
      {
        string const conditionKey = useSurfaceConditions ? "surface" : "reservoir";
        string const unitKey = useSurfaceConditions ? "s" : "r";

        real64 const currentTotalRate = connRate[iwelemRef];
        GEOS_LOG( GEOS_FMT( "{}: BHP (at the specified reference elevation): {} Pa",
                            wellControlsName, currentBHP ) );
        GEOS_LOG( GEOS_FMT( "{}: Total rate: {} kg/s; total {} volumetric rate: {} {}m3/s",
                            wellControlsName, currentTotalRate, conditionKey, currentTotalVolRate, unitKey ) );
        if( outputFile.is_open())
        {
          outputFile << "," << currentBHP;
          outputFile << "," << currentTotalRate << "," << currentTotalVolRate << std::endl;
          outputFile.close();
        }
      } );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseWell, string const &, Group * const )
}// namespace geos
