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

/**
 * @file surfaceGenerationKernels.hpp
 */

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/SolidBase.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "surfaceGenerationKernelsHelpers.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

namespace geos
{

namespace surfaceGenerationKernels
{

using namespace surfaceGenerationKernelsHelpers;

class NodalForceKernel
{

public:

  NodalForceKernel( ElementRegionManager const & elemManager,
                    constitutive::ConstitutiveManager const & constitutiveManager,
                    string const solidMaterialKey ):
    m_dNdX( elemManager.constructViewAccessor< array4d< real64 >, arrayView4d< real64 const > >( dataRepository::keys::dNdX ) ),
    m_detJ( elemManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( dataRepository::keys::detJ ) ),
    m_bulkModulus( elemManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "bulkModulus", constitutiveManager ) ),
    m_shearModulus( elemManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "shearModulus", constitutiveManager ) ),
    m_stress( elemManager.constructFullMaterialViewAccessor< array3d< real64, solid::STRESS_PERMUTATION >,
                                                             arrayView3d< real64 const, solid::STRESS_USD > >( constitutive::SolidBase::viewKeyStruct::stressString(),
                                                                                                               constitutiveManager ) )
  {
    m_solidMaterialFullIndex.resize( elemManager.numRegions() );
    elemManager.forElementRegionsComplete< CellElementRegion >( [&]( localIndex regionIndex,
                                                                     CellElementRegion const & region )
    {
      string const & solidMaterialName = region.getSubRegion( 0 ).getReference< string >( solidMaterialKey );
      constitutive::ConstitutiveBase const & solid = constitutiveManager.getConstitutiveRelation< constitutive::ConstitutiveBase >( solidMaterialName );
      m_solidMaterialFullIndex[regionIndex] = solid.getIndexInParent();
    } );
  }

  virtual void
  calculateSingleNodalForce( localIndex const er,
                             localIndex const esr,
                             localIndex const ei,
                             localIndex const targetNode,
                             real64 ( & force )[ 3 ] ) const
  {
    GEOS_MARK_FUNCTION;

    localIndex const numQuadraturePoints = m_detJ[er][esr].size( 1 );

    // Loop over quadrature points
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      real64 const quadratureStress[6] = LVARRAY_TENSOROPS_INIT_LOCAL_6 ( m_stress[er][esr][m_solidMaterialFullIndex[er]][ei][q] );
      real64 const dNdX[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( m_dNdX[er][esr][ei][q][targetNode] );
      computeNodalForce( quadratureStress, dNdX, m_detJ[er][esr][ei][q], force );
    }

    //wu40: the nodal force need to be weighted by Young's modulus and possion's ratio.
    scaleNodalForce( m_bulkModulus[er][esr][m_solidMaterialFullIndex[er]][ei], m_shearModulus[er][esr][m_solidMaterialFullIndex[er]][ei], force );
  }

protected:

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > const m_dNdX;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const m_detJ;

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const m_bulkModulus;

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const m_shearModulus;

  ElementRegionManager::MaterialViewAccessor< arrayView3d< real64 const, solid::STRESS_USD > > const m_stress;

  array1d< integer > m_solidMaterialFullIndex;
};


class PoroElasticNodalForceKernel : public NodalForceKernel
{

public:
  PoroElasticNodalForceKernel( ElementRegionManager const & elemManager,
                               constitutive::ConstitutiveManager const & constitutiveManager,
                               string const solidMaterialKey,
                               string const porosityModelKey ):
    NodalForceKernel( elemManager, constitutiveManager, solidMaterialKey ),
    m_pressure( elemManager.constructArrayViewAccessor< real64, 1 >( fields::flow::pressure::key() ) ),
    m_biotCoefficient( elemManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "biotCoefficient", constitutiveManager ) )
  {
    m_porosityMaterialFullIndex.resize( elemManager.numRegions() );
    elemManager.forElementRegionsComplete< CellElementRegion >( [&]( localIndex regionIndex,
                                                                     CellElementRegion const & region )
    {
      string const & porosityModelName = region.getSubRegion( 0 ).getReference< string >( porosityModelKey );
      constitutive::ConstitutiveBase const & porosityModel = constitutiveManager.getConstitutiveRelation< constitutive::ConstitutiveBase >( porosityModelName );
      m_porosityMaterialFullIndex[regionIndex] = porosityModel.getIndexInParent();
    } );

  }

  void
  calculateSingleNodalForce( localIndex const er,
                             localIndex const esr,
                             localIndex const ei,
                             localIndex const targetNode,
                             real64 ( & force )[ 3 ] ) const override

  {
    GEOS_MARK_FUNCTION;

    localIndex const numQuadraturePoints = m_detJ[er][esr].size( 1 );

    // Loop over quadrature points
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      real64 totalStress[6] = LVARRAY_TENSOROPS_INIT_LOCAL_6 ( m_stress[er][esr][m_solidMaterialFullIndex[er]][ei][q] );
      real64 const dNdX[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( m_dNdX[er][esr][ei][q][targetNode] );
      /// TODO: make it work for the thermal case as well
      LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -m_biotCoefficient[er][esr][m_porosityMaterialFullIndex[er]][ei] * m_pressure[er][esr][ei] );

      computeNodalForce( totalStress, dNdX, m_detJ[er][esr][ei][q], force );
    }

    //wu40: the nodal force need to be weighted by Young's modulus and possion's ratio.
    scaleNodalForce( m_bulkModulus[er][esr][m_solidMaterialFullIndex[er]][ei], m_shearModulus[er][esr][m_solidMaterialFullIndex[er]][ei], force );
  }

private:

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const m_pressure;

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const m_biotCoefficient;

  array1d< integer > m_porosityMaterialFullIndex;

};

template< typename LAMBDA >
void kernelSelector( ElementRegionManager const & elemManager,
                     constitutive::ConstitutiveManager const & constitutiveManager,
                     string const solidMaterialKey,
                     integer const isPoroelastic,
                     LAMBDA && lambda )
{
  if( isPoroelastic == 0 )
  {
    lambda( NodalForceKernel( elemManager, constitutiveManager, solidMaterialKey ) );
  }
  else
  {
    string const porosityModelKey = constitutive::CoupledSolidBase::viewKeyStruct::porosityModelNameString();
    lambda( PoroElasticNodalForceKernel( elemManager, constitutiveManager, solidMaterialKey, porosityModelKey )  );
  }

}

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos
