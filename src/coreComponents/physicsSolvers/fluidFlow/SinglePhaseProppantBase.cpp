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
 * @file SinglePhaseProppantBase.cpp
 */


#include "SinglePhaseProppantBase.hpp"

#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidSelector.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/ProppantSolid.hpp"
#include "constitutive/solid/porosity/ProppantPorosity.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransportFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseProppantBaseKernels.hpp"

namespace geos
{

using namespace constitutive;

using constitutive::SlurryFluidBase;

template< typename POROUSWRAPPER_TYPE >
void execute3( POROUSWRAPPER_TYPE porousWrapper,
               SurfaceElementSubRegion & subRegion,
               arrayView1d< real64 const > const & oldHydraulicAperture,
               arrayView1d< real64 const > const & newHydraulicAperture,
               arrayView1d< real64 const > const & proppantPackVolumeFraction )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFromApertureAndProppantVolumeFraction( k, q,
                                                                      oldHydraulicAperture[k],
                                                                      newHydraulicAperture[k],
                                                                      proppantPackVolumeFraction[k] );
    }
  } );
}

SinglePhaseProppantBase::SinglePhaseProppantBase( const string & name,
                                                  Group * const parent ):
  SinglePhaseBase( name, parent )
{}

SinglePhaseProppantBase::~SinglePhaseProppantBase()
{}

void SinglePhaseProppantBase::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  string & fluidMaterialName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidMaterialName = PhysicsSolverBase::getConstitutiveName< SlurryFluidBase >( subRegion );
  GEOS_ERROR_IF( fluidMaterialName.empty(), GEOS_FMT( "{}: Fluid model not found on subregion {}",
                                                      getDataContext(), subRegion.getName() ) );
}

void SinglePhaseProppantBase::validateConstitutiveModels( DomainPartition & domain ) const
{
  // Validate fluid models in regions
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      fluidName = getConstitutiveName< SlurryFluidBase >( subRegion );
      GEOS_THROW_IF( fluidName.empty(),
                     GEOS_FMT( "{}: Fluid model not found on subregion {}",
                               getDataContext(), subRegion.getName() ),
                     InputError );
    } );
  } );
}

SinglePhaseBase::FluidPropViews SinglePhaseProppantBase::getFluidProperties( constitutive::ConstitutiveBase const & fluid ) const
{
  SlurryFluidBase const & slurryFluid = dynamicCast< SlurryFluidBase const & >( fluid );
  return { slurryFluid.density(),
           slurryFluid.dDensity_dPressure(),
           slurryFluid.viscosity(),
           slurryFluid.dViscosity_dPressure(),
           slurryFluid.getField< fields::singlefluid::density >().getDefaultValue(),
           slurryFluid.getField< fields::singlefluid::viscosity >().getDefaultValue() };
}

void SinglePhaseProppantBase::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const proppantConcentration = dataGroup.getField< fields::proppant::proppantConcentration >();
  arrayView2d< real64 const > const componentConcentration = dataGroup.getField< fields::proppant::componentConcentration >();
  arrayView2d< real64 const > const cellBasedFlux = dataGroup.getField< fields::proppant::cellBasedFlux >();
  arrayView1d< integer const > const isProppantBoundaryElement = dataGroup.getField< fields::proppant::isProppantBoundary >();

  string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  SlurryFluidBase & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, fluidName );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    singlePhaseProppantBaseKernels::FluidUpdateKernel::launch( fluidWrapper,
                                                               pres,
                                                               proppantConcentration,
                                                               componentConcentration,
                                                               cellBasedFlux,
                                                               isProppantBoundaryElement );
  } );
}


void SinglePhaseProppantBase::updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const proppantPackVolumeFraction = subRegion.getField< fields::proppant::proppantPackVolumeFraction >();

  arrayView1d< real64 const > const newHydraulicAperture = subRegion.getField< fields::flow::hydraulicAperture >();
  arrayView1d< real64 const > const oldHydraulicAperture = subRegion.getField< fields::flow::aperture0 >();

  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( solidName );

  constitutive::ConstitutivePassThru< ProppantSolid< ProppantPorosity, ProppantPermeability > >::execute( porousSolid, [=, &subRegion] ( auto & castedProppantSolid )
  {
    typename TYPEOFREF( castedProppantSolid ) ::KernelWrapper porousWrapper = castedProppantSolid.createKernelUpdates();

    execute3( porousWrapper, subRegion, newHydraulicAperture, oldHydraulicAperture, proppantPackVolumeFraction );

  } );

}
}
