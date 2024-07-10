/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PhaseFieldPoromechanicsSolver.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDPOROMECHANICSSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDPOROMECHANICSSOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/simplePDE/PhaseFieldDamageFEM.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geos
{

class PhaseFieldPoromechanicsSolver : public CoupledSolver< SinglePhasePoromechanics< SinglePhaseBase >, PhaseFieldDamageFEM >
{
public:

  using Base = CoupledSolver< SinglePhasePoromechanics< SinglePhaseBase >, PhaseFieldDamageFEM >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  PhaseFieldPoromechanicsSolver( const string & name,
                                 Group * const parent );

  ~PhaseFieldPoromechanicsSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "PhaseFieldPoromechanics";
  }

  string getCatalogName() const override { return catalogName(); }

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "PhaseFieldPoromechanics"; }

  enum class SolverType : integer
  {
    Poromechanics = 0,
    Damage = 1
  };

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override final;

  virtual void postInputInitialization() override final;

  /**
   * @brief accessor for the pointer to the poromechanics solver
   * @return a pointer to the poromechanics solver
   */
  SinglePhasePoromechanics< SinglePhaseBase > * poromechancisSolver() const
  {
    return std::get< toUnderlying( SolverType::Poromechanics ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  PhaseFieldDamageFEM * damageSolver() const
  {
    return std::get< toUnderlying( SolverType::Damage ) >( m_solvers );
  }

  virtual void mapSolutionBetweenSolvers( DomainPartition & Domain, integer const idx ) override final;

  void mapDamageAndGradientToQuadrature( DomainPartition & domain );

  void applyDamageOnTractionBC( DomainPartition & domain );

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override final {}

};

template< typename FE_TYPE >
struct DamageAndDamageGradientInterpolationKernel
{
  DamageAndDamageGradientInterpolationKernel( CellElementSubRegion const & subRegion ):
    m_numElems( subRegion.size() )
  {}

  void interpolateDamageAndGradient( arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes,
                                     arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const xNodes,
                                     arrayView1d< real64 const > const nodalDamage,
                                     arrayView2d< real64 > damageFieldOnMaterial,
                                     arrayView3d< real64 > damageGradOnMaterial )
  {
    forAll< parallelDevicePolicy<> >( m_numElems, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
      constexpr localIndex n_q_points = FE_TYPE::numQuadraturePoints;

      real64 xLocal[ numNodesPerElement ][ 3 ];
      real64 nodalDamageLocal[ numNodesPerElement ];

      for( localIndex a = 0; a < numNodesPerElement; ++a )
      {
        localIndex const localNodeIndex = elemToNodes( k, a );

        for( int dim=0; dim < 3; ++dim )
        {
          xLocal[a][dim] = xNodes[ localNodeIndex ][dim];
        }

        nodalDamageLocal[ a ] = nodalDamage[ localNodeIndex ];
      }

      for( localIndex q = 0; q < n_q_points; ++q )
      {
        real64 N[ numNodesPerElement ];
        FE_TYPE::calcN( q, N );

        real64 dNdX[ numNodesPerElement ][ 3 ];
        real64 const detJ = FE_TYPE::calcGradN( q, xLocal, dNdX );

        GEOS_UNUSED_VAR( detJ );

        real64 qDamage = 0.0;
        real64 qDamageGrad[3] = {0, 0, 0};
        FE_TYPE::valueAndGradient( N, dNdX, nodalDamageLocal, qDamage, qDamageGrad );

        damageFieldOnMaterial( k, q ) = qDamage;

        for( int dim=0; dim < 3; ++dim )
        {
          damageGradOnMaterial[k][q][dim] = qDamageGrad[dim];
        }
      }

    } );
  }

  localIndex m_numElems;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_PhaseFieldPoromechanicsSOLVER_HPP_ */
