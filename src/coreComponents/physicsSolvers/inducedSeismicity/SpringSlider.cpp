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
 * @file SpringSlider.cpp
 */

#include "SpringSlider.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mesh/DomainPartition.hpp"
#include "kernels/RateAndStateKernels.hpp"
#include "rateAndStateFields.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

template< typename RSSOLVER_TYPE >
SpringSlider< RSSOLVER_TYPE >::SpringSlider( const string & name,
                                             Group * const parent ):
  RSSOLVER_TYPE( name, parent )
{}

template< typename RSSOLVER_TYPE >
SpringSlider< RSSOLVER_TYPE >::~SpringSlider()
{
  // TODO Auto-generated destructor stub
}

template< typename RSSOLVER_TYPE >
void SpringSlider< RSSOLVER_TYPE >::registerDataOnMesh( Group & meshBodies )
{
  RSSOLVER_TYPE::registerDataOnMesh( meshBodies );

  this->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      SurfaceElementSubRegion & subRegion )
    {
      // 3-component functions on fault
      string const labels3Comp[3] = { "normal", "tangent1", "tangent2" };
      subRegion.registerField< contact::dispJump >( this->getName() ).
      setDimLabels( 1, labels3Comp ).
      reference().template resizeDimension< 1 >( 3 );
      subRegion.registerField< contact::traction >( this->getName() ).
        setDimLabels( 1, labels3Comp ).
        reference().template resizeDimension< 1 >( 3 );
      subRegion.registerField< contact::deltaSlip >( this->getName() ).
        setDimLabels( 1, labels3Comp ).
        reference().template resizeDimension< 1 >( 3 );  

      subRegion.registerWrapper< string >( viewKeyStruct::frictionLawNameString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );

      string & frictionLawName = subRegion.getReference< string >( viewKeyStruct::frictionLawNameString() );
      frictionLawName =PhysicsSolverBase::getConstitutiveName< FrictionBase >( subRegion );
      GEOS_ERROR_IF( frictionLawName.empty(), GEOS_FMT( "{}: FrictionBase model not found on subregion {}",
                                                        this->getDataContext(), subRegion.getDataContext() ) );
    } );
  } );
}

template< typename RSSOLVER_TYPE >
real64 SpringSlider< RSSOLVER_TYPE >::updateStresses( real64 const dt,
                                     DomainPartition & domain ) const
{
  // Spring-slider shear traction computation
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {

      arrayView2d< real64 const > const deltaSlip = subRegion.getField< fields::contact::deltaSlip >();
      arrayView2d< real64 > const traction        = subRegion.getField< fields::contact::traction >();

      string const & fricitonLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
      RateAndStateFriction const & frictionLaw = this->template getConstitutiveModel< RateAndStateFriction >( subRegion, fricitonLawName );

      RateAndStateFriction::KernelWrapper frictionKernelWrapper = frictionLaw.createKernelUpdates();

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

template class SpringSlider< ImplicitQDRateAndState >;

namespace
{
typedef SpringSlider< ImplicitQDRateAndState > ImplicitSpringSlider;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ImplicitSpringSlider, string const &, dataRepository::Group * const )
}

} // namespace geos
