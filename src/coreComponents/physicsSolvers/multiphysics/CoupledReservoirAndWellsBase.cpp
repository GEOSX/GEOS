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
 * @file CoupledReservoirAndWellsBase.cpp
 *
 */

#include "CoupledReservoirAndWellsBase.hpp"

namespace geos
{

namespace coupledReservoirAndWellsInternal
{

void
addCouplingNumNonzeros( SolverBase const * const solver,
                        DomainPartition & domain,
                        DofManager & dofManager,
                        arrayView1d< localIndex > const & rowLengths,
                        integer const resNumDof,
                        integer const wellNumDof,
                        string const & resElemDofName,
                        string const & wellElemDofName )
{
  solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                        MeshLevel const & meshLevel,
                                                                        arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = meshLevel.getElemManager();

    string const wellDofKey = dofManager.getKey( wellElemDofName );
    string const resDofKey = dofManager.getKey( resElemDofName );

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resElemDofNumber =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );

    ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const & resElemGhostRank =
      elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

    globalIndex const rankOffset = dofManager.rankOffset();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const, WellElementSubRegion const & subRegion )
    {
      PerforationData const * const perforationData = subRegion.getPerforationData();

      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      // get the well degrees of freedom and ghosting info
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get the well element indices corresponding to each perforation
      arrayView1d< localIndex const > const & perfWellElemIndex =
        perforationData->getField< fields::perforation::wellElementIndex >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const & resElementRegion =
        perforationData->getField< fields::perforation::reservoirElementRegion >();
      arrayView1d< localIndex const > const & resElementSubRegion =
        perforationData->getField< fields::perforation::reservoirElementSubRegion >();
      arrayView1d< localIndex const > const & resElementIndex =
        perforationData->getField< fields::perforation::reservoirElementIndex >();

      // Loop over perforations and increase row lengths for reservoir and well elements accordingly
      forAll< serialPolicy >( perforationData->size(), [=] ( localIndex const iperf )
      {
        // get the reservoir (sub)region and element indices
        localIndex const er = resElementRegion[iperf];
        localIndex const esr = resElementSubRegion[iperf];
        localIndex const ei = resElementIndex[iperf];
        localIndex const iwelem = perfWellElemIndex[iperf];

        if( resElemGhostRank[er][esr][ei] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( resElemDofNumber[er][esr][ei] - rankOffset );
          GEOS_ASSERT_GE( localRow, 0 );
          GEOS_ASSERT_GE( rowLengths.size(), localRow + resNumDof );

          for( integer idof = 0; idof < resNumDof; ++idof )
          {
            rowLengths[localRow + idof] += wellNumDof;
          }
        }

        if( wellElemGhostRank[iwelem] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( wellElemDofNumber[iwelem] - rankOffset );
          GEOS_ASSERT_GE( localRow, 0 );
          GEOS_ASSERT_GE( rowLengths.size(), localRow + wellNumDof );

          for( integer idof = 0; idof < wellNumDof; ++idof )
          {
            rowLengths[localRow + idof] += resNumDof;
          }
        }
      } );
    } );
  } );
}

bool validateWellPerforations( SolverBase const * const reservoirSolver,
                               WellSolverBase const * const wellSolver,
                               DomainPartition const & domain )
{
  std::pair< string, string > badPerforation;

  arrayView1d< string const > const flowTargetRegionNames =
    reservoirSolver->getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );

  wellSolver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                            MeshLevel const & meshLevel,
                                                                            arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = meshLevel.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const, WellElementSubRegion const & subRegion )
    {
      PerforationData const * const perforationData = subRegion.getPerforationData();
      WellControls const & wellControls = wellSolver->getWellControls( subRegion );

      arrayView1d< localIndex const > const & resElementRegion =
        perforationData->getField< fields::perforation::reservoirElementRegion >();

      // Loop over perforations and check the reservoir region to which each perforation is connected to
      // If the name of the region is not in the list of targetted regions, then we have a "bad" connection.
      for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
      {
        localIndex const er = resElementRegion[iperf];
        string const regionName = elemManager.getRegion( er ).getName();
        if( std::find( flowTargetRegionNames.begin(), flowTargetRegionNames.end(), regionName ) == flowTargetRegionNames.end())
        {
          badPerforation = {wellControls.getName(), regionName};
          break;
        }
      }
    } );
  } );

  localIndex const hasBadPerforations = MpiWrapper::max( badPerforation.first.empty() ? 0 : 1 );

  GEOS_THROW_IF( !badPerforation.first.empty(),
                 GEOS_FMT( "{}: The well {} has a connection to the region {} which is not targeted by the solver",
                           wellSolver->getDataContext(), badPerforation.first, badPerforation.second ),
                 std::runtime_error );
  return hasBadPerforations == 0;
}

}

} /* namespace geos */
