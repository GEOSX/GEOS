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
 * @file SpringSliderTractionUpdate.cpp
 */

#include "SpringSliderTractionUpdate.hpp"
#include "physicsSolvers/inducedSeismicity/QuasiDynamicEQ.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "constitutive/contact/RateAndStateFriction.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

namespace inducedSeismicity
{

using namespace dataRepository;
using namespace fields;

void SpringSliderTractionUpdate::registerMissingDataOnMesh( SurfaceElementSubRegion & subRegion, string const & solverName ) const
{
  string const labels2Comp[2] = {"tangent1", "tangent2" };
  // 3-component functions on fault
  string const labels3Comp[3] = { "normal", "tangent1", "tangent2" };
  subRegion.registerField< contact::dispJump >( solverName ).
    setDimLabels( 1, labels3Comp ).
    reference().resizeDimension< 1 >( 3 );
  subRegion.registerField< contact::traction >( solverName ).
    setDimLabels( 1, labels3Comp ).
    reference().resizeDimension< 1 >( 3 );
  subRegion.registerField< contact::deltaSlip >( solverName ).
    setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );

  subRegion.registerWrapper< string >( QuasiDynamicEQ::viewKeyStruct::frictionLawNameString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );
  
  std::get< 0 >( m_solvers )->setFrictionLawName( subRegion );
}

real64 SpringSliderTractionUpdate::updateFaultTraction( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                        real64 const & dt,
                                                        int const GEOS_UNUSED_PARAM( cycleNumber ),
                                                        DomainPartition & domain ) const
{
  GEOS_UNUSED_VAR( domain );
  // Spring-slider shear traction computation
  std::get< 0 >( m_solvers )->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                                           MeshLevel & mesh,
                                                                                           arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {

      arrayView2d< real64 const > const deltaSlip = subRegion.getField< contact::deltaSlip >();
      arrayView2d< real64 > const traction        = subRegion.getField< fields::contact::traction >();

      constitutive::RateAndStateFriction const & frictionLaw = std::get< 0 >( m_solvers )->getFrictionLaw( subRegion );

      constitutive::RateAndStateFriction::KernelWrapper frictionKernelWrapper = frictionLaw.createKernelUpdates();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        SpringSliderParameters springSliderParameters = SpringSliderParameters( traction[k][0],
                                                                                frictionKernelWrapper.getACoefficient( k ),
                                                                                frictionKernelWrapper.getBCoefficient( k ),
                                                                                frictionKernelWrapper.getDcCoefficient( k ) );


        traction[k][1] = traction[k][1] + springSliderParameters.tauRate * dt
                         - springSliderParameters.springStiffness * deltaSlip[k][0];
        traction[k][2] = traction[k][2] + springSliderParameters.tauRate * dt
                         - springSliderParameters.springStiffness * deltaSlip[k][1];
      } );
    } );
  } );

  return dt;
}

}

}
