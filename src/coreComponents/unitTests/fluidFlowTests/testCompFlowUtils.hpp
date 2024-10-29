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

#ifndef GEOS_TESTCOMPFLOWUTILS_HPP
#define GEOS_TESTCOMPFLOWUTILS_HPP

#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "mesh/MeshManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "testFlowUtils.hpp"

namespace geos
{

namespace testing
{

using namespace geos::constitutive;
using namespace geos::constitutive::multifluid;

void fillNumericalJacobian( arrayView1d< real64 const > const & residual,
                            arrayView1d< real64 const > const & residualOrig,
                            globalIndex const dofIndex,
                            real64 const eps,
                            CRSMatrixView< real64, globalIndex const > const & jacobian )
{
  forAll< parallelDevicePolicy<> >( residual.size(), [=] GEOS_HOST_DEVICE ( localIndex const row )
  {
    real64 const dRdX = ( residual[row] - residualOrig[row] ) / eps;
    if( fabs( dRdX ) > 0.0 )
    {
      jacobian.addToRow< parallelDeviceAtomic >( row, &dofIndex, &dRdX, 1 );
    }
  } );
}

void setupProblemFromXML( ProblemManager & problemManager, char const * const xmlInput )
{
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.loadString( xmlInput );
  if( !xmlResult )
  {
    GEOS_LOG_RANK_0( "XML parsed with errors!" );
    GEOS_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOS_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  int mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );

  dataRepository::Group & commandLine =
    problemManager.getGroup< dataRepository::Group >( problemManager.groupKeys.commandLine );

  commandLine.registerWrapper< integer >( problemManager.viewKeys.xPartitionsOverride.key() ).
    setApplyDefaultValue( mpiSize );

  xmlWrapper::xmlNode xmlProblemNode = xmlDocument.getChild( dataRepository::keys::ProblemManager );
  problemManager.processInputFileRecursive( xmlDocument, xmlProblemNode );

  DomainPartition & domain = problemManager.getDomainPartition();

  constitutive::ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( constitutiveManager.getName().c_str());
  constitutiveManager.processInputFileRecursive( xmlDocument, topLevelNode );

  MeshManager & meshManager = problemManager.getGroup< MeshManager >( problemManager.groupKeys.meshManager );
  meshManager.generateMeshLevels( domain );

  ElementRegionManager & elementManager = domain.getMeshBody( 0 ).getBaseDiscretization().getElemManager();
  topLevelNode = xmlProblemNode.child( elementManager.getName().c_str());
  elementManager.processInputFileRecursive( xmlDocument, topLevelNode );

  problemManager.problemSetup();
  problemManager.applyInitialConditions();
}

void testCompositionNumericalDerivatives( CompositionalMultiphaseFVM & solver,
                                          DomainPartition & domain,
                                          real64 const perturbParameter,
                                          real64 const relTol )
{
  integer const numComp = solver.numFluidComponents();

  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                         [&]( string const,
                                              MeshLevel & mesh,
                                              arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions( regionNames,
                                               [&]( localIndex const,
                                                    ElementSubRegionBase & subRegion )
    {
      SCOPED_TRACE( subRegion.getParent().getParent().getName() + "/" + subRegion.getName() );

      string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      arrayView1d< string const > const & components = fluid.componentNames();

      arrayView2d< real64, compflow::USD_COMP > const compDens =
        subRegion.getField< fields::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens_n =
        subRegion.getField< fields::flow::globalCompDensity_n >();

      arrayView2d< real64, compflow::USD_COMP > const compFrac =
        subRegion.getField< fields::flow::globalCompFraction >();

      arrayView3d< real64, compflow::USD_COMP_DC > const dCompFrac_dCompDens =
        subRegion.getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >();

      // reset the solver state to zero out variable updates
      solver.resetStateToBeginningOfStep( domain );

      // make a copy of unperturbed values of component fractions
      array2d< real64, compflow::LAYOUT_COMP > compFracOrig( subRegion.size(), numComp );
      compFracOrig.setValues< serialPolicy >( compFrac );

      // update component density and check derivatives
      for( integer jc = 0; jc < numComp; ++jc )
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb a single component density in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
          compDens[ei][jc] += dRho;
        } );

        // recompute global component fractions
        solver.updateGlobalComponentFraction( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          auto dZ_dRho = invertLayout( dCompFrac_dCompDens[ei].toSliceConst(), numComp, numComp );
          string var = "compDens[" + components[jc] + "]";

          real64 const delta = compDens[ei][jc] - compDens_n[ei][jc];
          checkDerivative( compFrac[ei].toSliceConst(),
                           compFracOrig[ei].toSliceConst(),
                           dZ_dRho[jc].toSliceConst(),
                           delta,
                           relTol,
                           "compFrac",
                           var,
                           components );
        } );
      }
    } );
  } );
}

void testPhaseVolumeFractionNumericalDerivatives( CompositionalMultiphaseFVM & solver,
                                                  DomainPartition & domain,
                                                  bool const isThermal,
                                                  real64 const perturbParameter,
                                                  real64 const relTol )
{
  using Deriv = multifluid::DerivativeOffset;

  integer const numComp = solver.numFluidComponents();
  integer const numPhase = solver.numFluidPhases();

  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                         [&]( string const,
                                              MeshLevel & mesh,
                                              arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions( regionNames,
                                               [&]( localIndex const,
                                                    ElementSubRegionBase & subRegion )
    {
      SCOPED_TRACE( subRegion.getParent().getParent().getName() + "/" + subRegion.getName() );

      string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseFVM::viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      arrayView1d< string const > const & components = fluid.componentNames();
      arrayView1d< string const > const & phases = fluid.phaseNames();

      arrayView1d< real64 > const pres =
        subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const pres_n =
        subRegion.getField< fields::flow::pressure_n >();

      arrayView2d< real64, compflow::USD_COMP > const compDens =
        subRegion.getField< fields::flow::globalCompDensity >();
      arrayView2d< real64, compflow::USD_COMP > const compDens_n =
        subRegion.getField< fields::flow::globalCompDensity_n >();

      arrayView2d< real64, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::flow::phaseVolumeFraction >();

      arrayView3d< real64, compflow::USD_PHASE_DC > const dPhaseVolFrac =
        subRegion.getField< fields::flow::dPhaseVolumeFraction >();

      // reset the solver state to zero out variable updates
      solver.resetStateToBeginningOfStep( domain );

      // make a copy of unperturbed values of component fractions
      array2d< real64, compflow::LAYOUT_PHASE > phaseVolFracOrig( subRegion.size(), numPhase );
      phaseVolFracOrig.setValues< serialPolicy >( phaseVolFrac );

      // Step 1: update pressure and check derivatives

      {
        // perturb pressure in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseVolFracOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );
          auto dS = invertLayout( dPhaseVolFrac[ei].toSliceConst(), numPhase, numComp + 2 );

          real64 const delta = pres[ei] - pres_n[ei];
          checkDerivative( phaseVolFrac[ei].toSliceConst(),
                           phaseVolFracOrig[ei].toSliceConst(),
                           dS[Deriv::dP].toSliceConst(),
                           delta,
                           relTol,
                           "phaseVolFrac",
                           "Pres",
                           phases );
        } );
      }

      // Step 2: update component density and check derivatives

      for( integer jc = 0; jc < numComp; ++jc )
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb a single component density in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
          compDens[ei][jc] += dRho;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseVolFracOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          auto dS = invertLayout( dPhaseVolFrac[ei].toSliceConst(), numPhase, numComp + 2 );
          string var = "compDens[" + components[jc] + "]";

          real64 const delta = compDens[ei][jc] - compDens_n[ei][jc];
          checkDerivative( phaseVolFrac[ei].toSliceConst(),
                           phaseVolFracOrig[ei].toSliceConst(),
                           dS[Deriv::dC+jc].toSliceConst(),
                           delta,
                           relTol,
                           "phaseVolFrac",
                           var,
                           phases );
        } );
      }

      // Step 3: update temperature and check derivatives

      if( isThermal )
      {
        arrayView1d< real64 > const temp =
          subRegion.getField< fields::flow::temperature >();
        arrayView1d< real64 > const temp_n =
          subRegion.getField< fields::flow::temperature_n >();

        // reset the solver state to zero out variable updates
        solver.resetStateToBeginningOfStep( domain );

        // perturb temperature in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dT = perturbParameter * ( temp[ei] + perturbParameter );
          temp[ei] += dT;
        } );

        // recompute all fluid properties
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseVolFracOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          auto dS = invertLayout( dPhaseVolFrac[ei].toSliceConst(), numPhase, numComp + 2 );

          real64 const delta = temp[ei] - temp_n[ei];
          checkDerivative( phaseVolFrac[ei].toSliceConst(),
                           phaseVolFracOrig[ei].toSliceConst(),
                           dS[Deriv::dT].toSliceConst(),
                           delta,
                           relTol,
                           "phaseVolFrac",
                           "Temp",
                           phases );
        } );
      }

    } );
  } );
}

void testPhaseMobilityNumericalDerivatives( CompositionalMultiphaseFVM & solver,
                                            DomainPartition & domain,
                                            bool const isThermal,
                                            real64 const perturbParameter,
                                            real64 const relTol )
{
  using Deriv = multifluid::DerivativeOffset;

  integer const numComp = solver.numFluidComponents();
  integer const numPhase = solver.numFluidPhases();

  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                         [&]( string const,
                                              MeshLevel & mesh,
                                              arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions( regionNames,
                                               [&]( localIndex const,
                                                    ElementSubRegionBase & subRegion )
    {
      SCOPED_TRACE( subRegion.getParent().getParent().getName() + "/" + subRegion.getName() );

      string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseFVM::viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      arrayView1d< string const > const & components = fluid.componentNames();
      arrayView1d< string const > const & phases = fluid.phaseNames();

      arrayView1d< real64 > const pres =
        subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const pres_n =
        subRegion.getField< fields::flow::pressure_n >();

      arrayView2d< real64, compflow::USD_COMP > const compDens =
        subRegion.getField< fields::flow::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens_n =
        subRegion.getField< fields::flow::globalCompDensity_n >();

      arrayView2d< real64, compflow::USD_PHASE > const phaseMob =
        subRegion.getField< fields::flow::phaseMobility >();

      arrayView3d< real64, compflow::USD_PHASE_DC > const dPhaseMob =
        subRegion.getField< fields::flow::dPhaseMobility >();

      // reset the solver state to zero out variable updates
      solver.resetStateToBeginningOfStep( domain );

      // make a copy of unperturbed values of component fractions
      array2d< real64, compflow::LAYOUT_PHASE > phaseMobOrig( subRegion.size(), numPhase );
      phaseMobOrig.setValues< serialPolicy >( phaseMob );

      // Step 1: update pressure and check derivatives

      {
        // perturb pressure in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseMobOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          auto dMob = invertLayout( dPhaseMob[ei].toSliceConst(), numPhase, numComp + 2 );

          real64 const delta = pres[ei] - pres_n[ei];
          checkDerivative( phaseMob[ei].toSliceConst(),
                           phaseMobOrig[ei].toSliceConst(),
                           dMob[Deriv::dP].toSliceConst(),
                           delta,
                           relTol,
                           "phaseMob",
                           "Pres",
                           phases );
        } );
      }

      // Step 2: update component density and check derivatives

      for( integer jc = 0; jc < numComp; ++jc )
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb a single component density in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
          compDens[ei][jc] += dRho;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseMobOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          auto dMob = invertLayout( dPhaseMob[ei].toSliceConst(), numPhase, numComp + 2 );
          string var = "compDens[" + components[jc] + "]";

          real64 const delta = compDens[ei][jc] - compDens_n[ei][jc];
          checkDerivative( phaseMob[ei].toSliceConst(),
                           phaseMobOrig[ei].toSliceConst(),
                           dMob[Deriv::dC+jc].toSliceConst(),
                           delta,
                           relTol,
                           "phaseMob",
                           var,
                           phases );
        } );
      }

      // Step 3: update temperature and check derivatives

      if( isThermal )
      {
        arrayView1d< real64 > const temp =
          subRegion.getField< fields::flow::temperature >();
        arrayView1d< real64 const > const temp_n =
          subRegion.getField< fields::flow::temperature_n >();

        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb temperature in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dT = perturbParameter * ( temp[ei] + perturbParameter );
          temp[ei] += dT;
        } );

        // recompute component fractions
        solver.updateFluidState( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=, &phaseMobOrig] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          auto dMob = invertLayout( dPhaseMob[ei].toSliceConst(), numPhase, numComp + 2 );

          real64 const delta = temp[ei] - temp_n[ei];
          checkDerivative( phaseMob[ei].toSliceConst(),
                           phaseMobOrig[ei].toSliceConst(),
                           dMob[Deriv::dT].toSliceConst(),
                           delta,
                           relTol,
                           "phaseMob",
                           "Temp",
                           phases );
        } );
      }

    } );
  } );
}

template< typename COMPOSITIONAL_SOLVER, typename LAMBDA >
void fillCellCenteredNumericalJacobian( COMPOSITIONAL_SOLVER & solver,
                                        DomainPartition & domain,
                                        bool const isThermal,
                                        real64 const perturbParameter,
                                        arrayView1d< real64 > residual,
                                        arrayView1d< real64 > residualOrig,
                                        CRSMatrixView< real64, globalIndex > jacobian,
                                        CRSMatrixView< real64, globalIndex > jacobianFD,
                                        LAMBDA assembleFunction )
{
  integer const numComp = solver.numFluidComponents();

  DofManager const & dofManager = solver.getDofManager();
  string const elemDofKey = dofManager.getKey( COMPOSITIONAL_SOLVER::viewKeyStruct::elemDofFieldString() );

  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                       MeshLevel & mesh,
                                                                       arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const & elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );

      arrayView1d< real64 > const pres =
        subRegion.getField< fields::flow::pressure >();
      arrayView2d< real64, compflow::USD_COMP > const compDens =
        subRegion.getField< fields::flow::globalCompDensity >();

      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        if( elemGhostRank[ei] >= 0 )
        {
          continue;
        }

        // Warning! To compute a correct totalDensity, we have to reset compDens
        solver.resetStateToBeginningOfStep( domain );

        real64 totalDensity = 0.0;
        compDens.move( hostMemorySpace, true ); // to get the correct compDens after reset
        for( integer ic = 0; ic < numComp; ++ic )
        {
          totalDensity += compDens[ei][ic];
        }

        // Step 1: compute numerical derivatives wrt pressure

        {
          pres.move( hostMemorySpace, true ); // to get the correct pres after reset
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
#if defined(GEOS_USE_CUDA)
          pres.move( parallelDeviceMemorySpace, false );
#endif


          mesh.getElemManager().forElementSubRegions( regionNames,
                                                      [&]( localIndex const,
                                                           ElementSubRegionBase & subRegion2 )
          {
            solver.updateFluidState( subRegion2 );
          } );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 elemDofNumber[ei],
                                 dP,
                                 jacobianFD.toViewConstSizes() );
        }

        // Step 2: compute numerical derivatives wrt component densities

        for( integer jc = 0; jc < numComp; ++jc )
        {
          solver.resetStateToBeginningOfStep( domain );

          compDens.move( hostMemorySpace, true ); // to get the correct compDens after reset
          real64 const dRho = perturbParameter * totalDensity;
          compDens[ei][jc] += dRho;
#if defined(GEOS_USE_CUDA)
          compDens.move( parallelDeviceMemorySpace, false );
#endif


          mesh.getElemManager().forElementSubRegions( regionNames,
                                                      [&]( localIndex const,
                                                           ElementSubRegionBase & subRegion2 )
          {
            solver.updateFluidState( subRegion2 );
          } );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 elemDofNumber[ei] + jc + 1,
                                 dRho,
                                 jacobianFD.toViewConstSizes() );
        }

        // Step 3: compute numerical derivatives wrt temperature

        if( isThermal )
        {
          arrayView1d< real64 > const temp =
            subRegion.getField< fields::flow::temperature >();
          temp.move( hostMemorySpace, false );

          solver.resetStateToBeginningOfStep( domain );

          temp.move( hostMemorySpace, true );
          real64 const dT = perturbParameter * ( temp[ei] + perturbParameter );
          temp[ei] += dT;
#if defined(GEOS_USE_CUDA)
          temp.move( parallelDeviceMemorySpace, false );
#endif

          // here, we make sure that rock internal energy is updated
          // in other words, a call to updateFluidState would not work
          solver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 elemDofNumber[ei] + numComp + 1,
                                 dT,
                                 jacobianFD.toViewConstSizes() );
        }

      }
    } );
  } );
}


} // namespace testing

} // namespace geos

#endif //GEOS_TESTCOMPFLOWUTILS_HPP
