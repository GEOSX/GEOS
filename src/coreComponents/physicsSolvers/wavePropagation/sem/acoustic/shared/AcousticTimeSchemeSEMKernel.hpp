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
 * @file AcousticTimeSchemeSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICTIMESCHEMESEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICTIMESCHEMESEMKERNEL_HPP_

namespace geos
{

struct AcousticTimeSchemeSEM
{

  using EXEC_POLICY = parallelDevicePolicy< >;


  /**
   * @brief  Apply second order Leap-Frog time scheme for isotropic case without PML
   * @param[in] dt time-step
   * @param[out] p_np1 pressure array at time n+1 (updated here)
   * @param[in] p_n pressure array at time n
   * @param[in] p_nm1 pressure array at time n-1
   * @param[in] mass the mass matrix
   * @param[in] stiffnessVector array containing the product of the stiffness matrix R and the pressure at time n
   * @param[in] damping the damping matrix
   * @param[in] rhs the right-hand-side
   * @param[in] freeSurfaceNodeIndicator array which contains indicators to tell if we are on a free-surface boundary or not
   * @param[in] solverTargetNodesSet the targetted nodeset (useful in particular when we do elasto-acoustic simulation )
   */
  static void LeapFrogWithoutPML( real64 const dt,
                                  arrayView1d< real32 > const p_np1,
                                  arrayView1d< real32 > const p_n,
                                  arrayView1d< real32 > const p_nm1,
                                  arrayView1d< real32 const > const mass,
                                  arrayView1d< real32 > const stiffnessVector,
                                  arrayView1d< real32 const > const damping,
                                  arrayView1d< real32 > const rhs,
                                  arrayView1d< localIndex const > const freeSurfaceNodeIndicator,
                                  SortedArrayView< localIndex const > const solverTargetNodesSet )
  {
    real64 const dt2 = pow( dt, 2 );
    forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = solverTargetNodesSet[n];
      if( freeSurfaceNodeIndicator[a] != 1 )
      {
        p_np1[a] = p_n[a];
        p_np1[a] *= 2.0 * mass[a];
        p_np1[a] -= (mass[a] - 0.5 * dt * damping[a]) * p_nm1[a];
        p_np1[a] += dt2 * (rhs[a] - stiffnessVector[a]);
        p_np1[a] /= mass[a] + 0.5 * dt * damping[a];
      }
    } );

  };

  /**
   * @brief  Apply second order Leap-Frog time scheme for VTI case without PML
   * @param[in] size The number of nodes in the nodeManager
   * @param[in] dt time-step
   * @param[out] p_np1 pressure array at time n+1 (updated here)
   * @param[in] p_n pressure array at time n
   * @param[in] p_nm1 pressure array at time n-1
   * @param[out] q_np1 auxiliary pressure array at time n+1 (updated here)
   * @param[in] q_n auxiliary pressure array at time n
   * @param[in] q_nm1 auxiliary pressure array at time n-1
   * @param[in] mass the mass matrix
   * @param[in] stiffnessVector_p array containing the product of the stiffness matrix R and the pressure at time n
   * @param[in] stiffnessVector_q array containing the product of the stiffness matrix R and the auxiliary pressure at time n
   * @param[in] damping_p the damping matrix
   * @param[in] damping_pq the damping matrix
   * @param[in] damping_q the damping matrix
   * @param[in] damping_qp the damping matrix
   * @param[in] rhs the right-hand-side
   * @param[in] freeSurfaceNodeIndicator array which contains indicators to tell if we are on a free-surface boundary or not
   * @param[in] lateralSurfaceNodeIndicator array which contains indicators to tell if we are on a lateral boundary or not
   * @param[in] bottomSurfaceNodeIndicator array which contains indicators to telle if we are on the bottom boundary or not
   */
  static void LeapFrogforVTI( localIndex const size,
                              real64 const dt,
                              arrayView1d< real32 > const p_np1,
                              arrayView1d< real32 > const p_n,
                              arrayView1d< real32 > const p_nm1,
                              arrayView1d< real32 > const q_np1,
                              arrayView1d< real32 > const q_n,
                              arrayView1d< real32 > const q_nm1,
                              arrayView1d< real32 const > const mass,
                              arrayView1d< real32 > const stiffnessVector_p,
                              arrayView1d< real32 > const stiffnessVector_q,
                              arrayView1d< real32 const > const damping_p,
                              arrayView1d< real32 const > const damping_pq,
                              arrayView1d< real32 const > const damping_q,
                              arrayView1d< real32 const > const damping_qp,
                              arrayView1d< real32 > const rhs,
                              arrayView1d< localIndex const > const freeSurfaceNodeIndicator,
                              arrayView1d< localIndex const > const lateralSurfaceNodeIndicator,
                              arrayView1d< localIndex const > const bottomSurfaceNodeIndicator )

  {
    real64 const dt2 = pow( dt, 2 );
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      if( freeSurfaceNodeIndicator[a] != 1 )
      {
        p_np1[a] = 2.0*mass[a]*p_n[a]/dt2;
        p_np1[a] -= mass[a]*p_nm1[a]/dt2;
        p_np1[a] += stiffnessVector_p[a];
        p_np1[a] += rhs[a];

        q_np1[a] = 2.0*mass[a]*q_n[a]/dt2;
        q_np1[a] -= mass[a]*q_nm1[a]/dt2;
        q_np1[a] += stiffnessVector_q[a];
        q_np1[a] += rhs[a];

        if( lateralSurfaceNodeIndicator[a] != 1 && bottomSurfaceNodeIndicator[a] != 1 )
        {
          // Interior node, no boundary terms
          p_np1[a] /= mass[a]/dt2;
          q_np1[a] /= mass[a]/dt2;
        }
        else
        {
          // Boundary node
          p_np1[a] += damping_p[a]*p_nm1[a]/dt/2;
          p_np1[a] += damping_pq[a]*q_nm1[a]/dt/2;

          q_np1[a] += damping_q[a]*q_nm1[a]/dt/2;
          q_np1[a] += damping_qp[a]*p_nm1[a]/dt/2;
          // Hand-made Inversion of 2x2 matrix
          real32 coef_pp = mass[a]/dt2;
          coef_pp += damping_p[a]/dt/2;
          real32 coef_pq = damping_pq[a]/dt/2;

          real32 coef_qq = mass[a]/dt2;
          coef_qq += damping_q[a]/2/dt;
          real32 coef_qp = damping_qp[a]/dt/2;

          real32 det_pq = 1/(coef_pp * coef_qq - coef_pq*coef_qp);

          real32 aux_p_np1 = p_np1[a];
          p_np1[a] = det_pq*(coef_qq*p_np1[a] - coef_pq*q_np1[a]);
          q_np1[a] = det_pq*(coef_pp*q_np1[a] - coef_qp*aux_p_np1);
        }
      }
    } );
  };



};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICTIMESCHEMESEMKERNEL_HPP_
