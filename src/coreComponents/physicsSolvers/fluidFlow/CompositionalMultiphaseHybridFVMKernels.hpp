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
 * @file CompositionalMultiphaseHybridFVMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"


namespace geos
{

namespace compositionalMultiphaseHybridFVMKernels
{

// struct to specify local and neighbor derivatives
struct Pos
{
  static constexpr integer LOCAL = 0;
  static constexpr integer NEIGHBOR = 1;
};

/******************************** UpwindingHelper ********************************/

struct UpwindingHelper
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief At a given one-sided face, compute the upwind viscous transport coefficient
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens the derivatives of the phase densities in the domain wrt pressure and component fractions (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the domain (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac the derivatives of the phase component fractions in the domain wrt pressure and component fractions
   *(non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] oneSidedVolFlux the total volumetric flux at this face
   * @param[out] upwPhaseViscCoef the upwind viscous transport coef at this face
   * @param[out] dUpwPhaseViscCoef_dPres the derivative of the upwind viscous transport coef wrt pressure at this face
   * @param[out] dUpwPhaseViscCoef_dCompDens the derivatives of the upwind viscous transport coef wrt component density at this face
   * @param[out] upwViscDofNumber the dof number of the upwind cell at this face
   */
  template< integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    upwindViscousCoefficient( localIndex const (&localIds)[ 3 ],
                              localIndex const (&neighborIds)[ 3 ],
                              ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens,
                              ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseDens,
                              ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                              ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                              ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                              ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                              real64 const & oneSidedVolFlux,
                              real64 ( &upwPhaseViscCoef )[ NP ][ NC ],
                              real64 ( &dUpwPhaseViscCoef_dPres )[ NP ][ NC ],
                              real64 ( &dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ],
                              globalIndex & upwViscDofNumber );

  /**
   * @brief At a given one-sided face, compute the upwind viscous transport coefficient
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] transGravCoef
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens the derivatives of the phase densities in the domain wrt pressure and component fractions (non-local)
   * @param[in] phaseMassDens the phase mass densities in the domain (non-local)
   * @param[in] dPhaseMassDens the derivatives of the phase mass densities in the domain wrt pressure and component fractions (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the domain (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac the derivatives of the phase component fractions in the domain wrt pressure and component fractions
   *(non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[inout] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[inout] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   * @param[inout] upwPhaseGravCoef the upwinded buoyancy transport coefficient at this face (ifaceLoc)
   * @param[inout] dUpwPhaseGravCoef_dPres the derivative of the upwinded buoyancy transport coefficient wrt pressure
   * @param[inout] dUpwPhaseGravCoef_dCompDens the derivative of the upwinded buoyancy transport coefficient wrt component density
   */
  template< integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    upwindBuoyancyCoefficient( localIndex const (&localIds)[ 3 ],
                               localIndex const (&neighborIds)[ 3 ],
                               real64 const & transGravCoef,
                               ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens,
                               ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseDens,
                               ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                               ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                               ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                               ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                               ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                               ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                               ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                               real64 ( &phaseGravTerm )[ NP ][ NP-1 ],
                               real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                               real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ],
                               real64 ( &upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ],
                               real64 ( &dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ],
                               real64 ( &dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] );


  /**
   * @brief At a given one-sided face, compute the gravCoef multiplied by the difference in phase densities
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] transGravCoef
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens the derivatives of the phase densities in the domain wrt pressure and component fractions (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[inout] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[inout] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   */
  template< integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    computePhaseGravTerm( localIndex const (&localIds)[ 3 ],
                          localIndex const (&neighborIds)[ 3 ],
                          real64 const & transGravCoef,
                          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                          real64 ( &phaseGravTerm )[ NP ][ NP-1 ],
                          real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                          real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] );

  /**
   * @brief At a given one-sided face, compute the upwinded total mobility
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the domain (non-local)
   * @param[in] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] totalMob the upwinded total mobility
   * @param[inout] dTotalMob_dPres the derivative of the upwinded total mobility wrt pressure
   * @param[inout] dTotalMob_dCompDens the derivative of the upwinded total mobility wrt component density
   */
  template< integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    computeUpwindedTotalMobility( localIndex const (&localIds)[ 3 ],
                                  localIndex const (&neighborIds)[ 3 ],
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                  real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                                  real64 & totalMob,
                                  real64 ( &dTotalMob_dPres )[ 2 ],
                                  real64 ( &dTotalMob_dCompDens )[ 2 ][ NC ] );

  /**
   * @brief Set the element indices used to evaluate the mobility ratios of the buoyancy term in hybrid upwinding
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] gravTerm the gravCoef multiplied by the difference in phase densities
   * @param[in] totalMob the upwinded total mobility
   * @param[in] eru region index of the upwind element
   * @param[in] esru subRegion index of the upwind element
   * @param[in] eiu element index of the upwind element
   * @param[in] posu position (local or neighbor) of the upwind element
   * @param[in] erd region index of the downwind element
   * @param[in] esrd subRegion index of the downwind element
   * @param[in] eid element index of the downwind element
   * @param[in] posd position (local or neighbor) of the downwind element
   */
  GEOS_HOST_DEVICE
  inline
  static
  void
  setIndicesForMobilityRatioUpwinding( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const & gravTerm,
                                       localIndex & eru, localIndex & esru, localIndex & eiu, localIndex & posu,
                                       localIndex & erd, localIndex & esrd, localIndex & eid, localIndex & posd );

  /**
   * @brief Set the element indices used to evaluate the total mobility of the buoyancy term in hybrid upwinding
   * @param[in] localIds triplet of indices for the local element
   * @param[in] neighborIds triplet of indices for the neighbor element
   * @param[in] gravTerm gravity term used to upwind
   * @param[out] totalMobIds for each phase, triplet of indices of the upwind element
   * @param[out] totalMobPos for each phase, flag specifying with the upwind element is local or neighbor
   */
  template< integer NP >
  GEOS_HOST_DEVICE
  static void
    setIndicesForTotalMobilityUpwinding( localIndex const (&localIds)[ 3 ],
                                         localIndex const (&neighborIds)[ 3 ],
                                         real64 const (&gravTerm)[ NP ][ NP-1 ],
                                         localIndex ( &totalMobIds )[ NP ][ 3 ],
                                         localIndex ( &totalMobPos )[ NP ] );

};


/******************************** AssemblerKernelHelper ********************************/

struct AssemblerKernelHelper
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief In a given element, compute the transmissibility-weighted pressure gradients in the cell
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] faceGravCoef the depth at the mesh facesb
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] phaseDens the phase densities in the element
   * @param[in] dPhaseDens the derivatives of the phase densities in the element wrt pressure and component fractions
   * @param[in] phaseMob the phase mobilities in the element
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the element
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the element wrt component density
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[out] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[out] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[out] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    applyGradient( arrayView1d< real64 const > const & facePres,
                   arrayView1d< real64 const > const & faceGravCoef,
                   arraySlice1d< localIndex const > const & elemToFaces,
                   real64 const & elemPres,
                   real64 const & elemGravCoef,
                   arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & elemPhaseMassDens,
                   arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const & dElemPhaseMassDens,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & elemPhaseMob,
                   arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dElemPhaseMob,
                   arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dElemCompFrac_dCompDens,
                   arraySlice2d< real64 const > const & transMatrix,
                   real64 ( &oneSidedVolFlux )[ NF ],
                   real64 ( &dOneSidedVolFlux_dPres )[ NF ],
                   real64 ( &dOneSidedVolFlux_dFacePres )[ NF ][ NF ],
                   real64 ( &dOneSidedVolFlux_dCompDens )[ NF ][ NC ] );

  /**
   * @brief In a given element, compute the flux divergence, i.e, sum the fluxes at this element's faces
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] rankOffset offset of this rank
   * @param[in] faceDofNumber the dof numbers of the face pressures in the domain
   * @param[in] elemRegionList map from face to element region index
   * @param[in] elemSubRegionList map from face to element subRegion index
   * @param[in] elemList map from face to element index
   * @param[in] regionFilter set of target regions of the solver
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] phaseDens phaseDens the phase densities in the domain
   * @param[in] dPhaseDens the derivatives of the phase densities wrt pressure and component fractions
   * @param[in] phaseMob the phase mobilities in the domain
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the element
   * @param[in] dCompFrac_dCompDens the derivatives of the component fractions wrt component density
   * @param[in] phaseCompFrac the phase component fractions in the domain
   * @param[in] dPhaseCompFrac the derivatives of the phase component fractions wrt pressure and component fractions
   * @param[in] elemDofNumber the dof numbers of the element in the domain
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[in] dt the time step
   * @param[inout] localMatrix the Jacobian matrix
   * @param[inout] localRhs the residual
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
  assembleFluxDivergence( localIndex const (&localIds)[ 3 ],
                          globalIndex const rankOffset,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          SortedArrayView< localIndex const > const & regionFilter,
                          arrayView1d< globalIndex const > const & faceDofNumber,
                          arrayView1d< real64 const > const & mimFaceGravCoef,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          real64 const & elemGravCoef,
                          integer const useTotalMassEquation,
                          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens,
                          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseDens,
                          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                          ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                          arraySlice2d< real64 const > const & transMatrixGrav,
                          real64 const (&oneSidedVolFlux)[ NF ],
                          real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs );

  /**
   * @brief In a given element, compute the viscous flux divergence, i.e, sum the viscous fluxes at this element's faces
   * @param[in] ifaceLoc the local index of the face
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face centered pressure
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[in] upwPhaseViscCoef the upwinded viscous transport coefficient at this face (ifaceLoc)
   * @param[in] dUpwPhaseViscCoef_dPres the derivative of the upwinded viscous transport coefficient wrt pressure
   * @param[in] dUpwPhaseViscCoef_dCompDens the derivative of the upwinded viscous transport coefficient wrt component density
   * @param[in] upwViscDofNumber degree of freedom number of the upwind element
   * @param[in] faceDofNumber degree of freedom number of the face
   * @param[in] dt the time step
   * @param[inout] divMassFluxes the divergence of the fluxes in the element
   * @param[inout] dDivMassFluxes_dElemVars the derivatives of the flux divergence wrt the element centered vars (pres and comp dens)
   * @param[inout] dDivMassFluxes_dFaceVars the derivatives of the flux divergence wrt the face centered vars
   * @param[inout] dofColIndicesElemVars degrees of freedom of the cells involved in the flux divergence
   * @param[inout] dofColIndicesFaceVars degrees of freedom of the faces involved in the flux divergence
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    assembleViscousFlux( localIndex const ifaceLoc,
                         real64 const (&oneSidedVolFlux)[ NF ],
                         real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                         real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                         real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                         real64 const (&upwPhaseViscCoef)[ NP ][ NC ],
                         real64 const (&dUpwPhaseViscCoef_dPres)[ NP ][ NC ],
                         real64 const (&dUpwPhaseViscCoef_dCompDens)[ NP ][ NC ][ NC ],
                         globalIndex const elemDofNumber,
                         globalIndex const neighborDofNumber,
                         globalIndex const upwViscDofNumber,
                         globalIndex const faceDofNumber,
                         real64 const & dt,
                         real64 ( &divMassFluxes )[ NC ],
                         real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ],
                         real64 ( &dDivMassFluxes_dFaceVars )[ NC ][ NF ],
                         globalIndex ( &dofColIndicesElemVars )[ (NC+1)*(NF+1) ],
                         globalIndex ( &dofColIndicesFaceVars )[ NF ] );

  /**
   * @brief In a given element, compute the buoyancy flux divergence, i.e, sum the buoyancy fluxes at this element's faces
   * @param[in] ifaceLoc the local index of the face
   * @param[in] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[in] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[in] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   * @param[in] upwPhaseGravCoef the upwinded buoyancy transport coefficient at this face (ifaceLoc)
   * @param[in] dUpwPhaseGravCoef_dPres the derivative of the upwinded buoyancy transport coefficient wrt pressure
   * @param[in] dUpwPhaseGravCoef_dCompDens the derivative of the upwinded buoyancy transport coefficient wrt component density
   * @param[in] neighborDofNumber the degree of freedom number of the neighbor element
   * @param[in] dt the time step
   * @param[inout] divMassFluxes the divergence of the fluxes in the element
   * @param[inout] dDivMassFluxes_dElemVars the derivatives of the flux divergence wrt the element centered vars (pres and comp dens)
   * @param[inout] dofColIndicesElemVars degrees of freedom of the cells involved in the flux divergence
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
    assembleBuoyancyFlux( localIndex const ifaceLoc,
                          real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                          real64 const (&dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ],
                          real64 const (&dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ],
                          real64 const (&upwPhaseGravCoef)[ NP ][ NP-1 ][ NC ],
                          real64 const (&dUpwPhaseGravCoef_dPres)[ NP ][ NP-1 ][ NC ][ 2 ],
                          real64 const (&dUpwPhaseGravCoef_dCompDens)[ NP ][ NP-1 ][ NC ][ 2 ][ NC ],
                          real64 const & dt,
                          real64 ( &divMassFluxes )[ NC ],
                          real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ] );

  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures in the subRegion
   * @param[in] faceGhostRank the ghost ranks of the face pressures in the subRegion
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] rankOffset the offset of this rank
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
  assembleFaceConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                           arrayView1d< integer const > const & faceGhostRank,
                           arraySlice1d< localIndex const > const & elemToFaces,
                           globalIndex const elemDofNumber,
                           globalIndex const rankOffset,
                           real64 const (&oneSidedVolFlux)[ NF ],
                           real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                           real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                           real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs );

};


/******************************** AssemblerKernel ********************************/

struct AssemblerKernel
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief In a given element, assemble the mass conservation equations and the contribution of this element to the face
   * constraints
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] faceGhostRank ghost rank of each face
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens the derivatives of the phase densities in the domain wrt pressure and component fractions (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the domain (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac the derivatives of the phase component fractions in the domain wrt pressure and component fractions
   *(non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] elemGhostRank the ghost rank of the element in which we assemble the fluxes
   * @param[in] rankOffset the offset of this rank
   * @param[in] lengthTolerance tolerance used in the transmissibility matrix computation
   * @param[in] dt time step size
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< integer NF, integer NC, integer NP >
  GEOS_HOST_DEVICE
  static void
  compute( localIndex const er, localIndex const esr, localIndex const ei,
           SortedArrayView< localIndex const > const & regionFilter,
           arrayView2d< localIndex const > const & elemRegionList,
           arrayView2d< localIndex const > const & elemSubRegionList,
           arrayView2d< localIndex const > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< integer const > const & faceGhostRank,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arrayView1d< real64 const > const & mimFaceGravCoef,
           arraySlice1d< localIndex const > const & elemToFaces,
           real64 const & elemPres,
           real64 const & elemGravCoef,
           const integer useTotalMassEquation,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseDens,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           arraySlice2d< real64 const > const & transMatrixGrav,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs );
};


/******************************** FluxKernel ********************************/

struct FluxKernel
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using CompFlowAccessors =
    StencilAccessors< fields::flow::phaseMobility,
                      fields::flow::dPhaseMobility,
                      fields::flow::dGlobalCompFraction_dGlobalCompDensity >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseMassDensity,
                              fields::multifluid::dPhaseMassDensity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;


  /**
   * @brief In a given subRegion, assemble the mass conservation equations and the contribution of the elements of this subRegion  to the
   * face constraints.
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion the subRegion in which we are going to assemble the fluxes
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] faceGhostRank  the ghost ranks of the face pressures
   * @param[in] facePres the pressures at the mesh faces at the beginning of the time step
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens the derivatives of the phase densities in the domain wrt pressure and component fractions (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob the derivatives of the phase mobilities in the domain (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac the derivatives of the phase component fractions in the domain wrt pressure and component fractions
   *(non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] rankOffset the offset of this rank
   * @param[in] lengthTolerance tolerance used in the transmissibility matrix computation
   * @param[in] dt time step size
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< integer NF, integer NC, integer NP, typename IP_TYPE >
  static void
  launch( localIndex er, localIndex esr,
          CellElementSubRegion const & subRegion,
          constitutive::PermeabilityBase const & permeabilityModel,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView1d< globalIndex const > const & faceDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & faceGravCoef,
          arrayView1d< real64 const > const & mimFaceGravCoef,
          arrayView1d< real64 const > const & transMultiplier,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseDens,
          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
          globalIndex const rankOffset,
          real64 const lengthTolerance,
          real64 const dt,
          integer const useTotalMassEquation,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};


/******************************** PhaseMobilityKernel ********************************/

/**
 * @class PhaseMobilityKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase mobilities
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseMobilityKernel : public isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  PhaseMobilityKernel( ObjectManagerBase & subRegion,
                       constitutive::MultiFluidBase const & fluid,
                       constitutive::RelativePermeabilityBase const & relperm )
    : Base(),
    m_phaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::flow::dPhaseVolumeFraction >() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseVisc( fluid.phaseViscosity() ),
    m_dPhaseVisc( fluid.dPhaseViscosity() ),
    m_phaseRelPerm( relperm.phaseRelPerm() ),
    m_dPhaseRelPerm_dPhaseVolFrac( relperm.dPhaseRelPerm_dPhaseVolFraction() ),
    m_phaseMob( subRegion.getField< fields::flow::phaseMobility >() ),
    m_dPhaseMob( subRegion.getField< fields::flow::dPhaseMobility >() )
  {}

  /**
   * @brief Compute the phase mobilities in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseMobilityKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && phaseMobilityKernelOp = NoOpFunc{} ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseVisc = m_phaseVisc[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseVisc = m_dPhaseVisc[ei][0];
    arraySlice1d< real64 const, constitutive::relperm::USD_RELPERM - 2 > const phaseRelPerm = m_phaseRelPerm[ei][0];
    arraySlice2d< real64 const, constitutive::relperm::USD_RELPERM_DS - 2 > const dPhaseRelPerm_dPhaseVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[ei][0];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dPhaseVolFrac[ei];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const phaseMob = m_phaseMob[ei];
    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseMob = m_dPhaseMob[ei];

    real64 dRelPerm_dC[numComp]{};
    real64 dVisc_dC[numComp]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      // compute the phase mobility only if the phase is present
      bool const phaseExists = (phaseVolFrac[ip] > 0);
      if( !phaseExists )
      {
        phaseMob[ip] = 0.;
        for( integer jc = 0; jc < numComp + 2; ++jc )
        {
          dPhaseMob[ip][jc] = 0.;
        }
        continue;
      }

      real64 const viscosity = phaseVisc[ip];
      real64 const dVisc_dP = dPhaseVisc[ip][Deriv::dP];
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseVisc[ip], dVisc_dC, Deriv::dC );

      real64 const relPerm = phaseRelPerm[ip];
      real64 dRelPerm_dP = 0.0;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dRelPerm_dC[ic] = 0.0;
      }

      for( integer jp = 0; jp < numPhase; ++jp )
      {
        real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
        dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac[jp][Deriv::dP];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac[jp][Deriv::dC+jc];
        }
      }

      real64 const mobility = relPerm / viscosity;

      phaseMob[ip] = mobility;
      dPhaseMob[ip][Deriv::dP] = dRelPerm_dP / viscosity
                                 - mobility * dVisc_dP / viscosity;

      // compositional derivatives
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseMob[ip][Deriv::dC+jc] = dRelPerm_dC[jc] / viscosity
                                      - mobility * dVisc_dC[jc] / viscosity;
      }

      // call the lambda in the phase loop to allow the reuse of the relperm, viscosity, and mobility
      // possible use: assemble the derivatives wrt temperature
      phaseMobilityKernelOp( ip, phaseMob[ip], dPhaseMob[ip] );
    }
  }


protected:

  // inputs

  /// Views on the phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dPhaseVolFrac;
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

  /// Views on the phase viscosities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseVisc;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseVisc;

  /// Views on the phase relative permeabilities
  arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > m_phaseRelPerm;
  arrayView4d< real64 const, constitutive::relperm::USD_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;

  // outputs

  /// Views on the phase mobilities
  arrayView2d< real64, compflow::USD_PHASE > m_phaseMob;
  arrayView3d< real64, compflow::USD_PHASE_DC > m_dPhaseMob;

};

/**
 * @class PhaseMobilityKernelFactory
 */
class PhaseMobilityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   constitutive::MultiFluidBase const & fluid,
                   constitutive::RelativePermeabilityBase const & relperm )
  {
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 2 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 3 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 1 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using CompFlowAccessors =
    StencilAccessors< fields::elementVolume >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::totalDensity_n >;
  using PorosityAccessors =
    StencilMaterialAccessors< constitutive::PorosityBase,
                              fields::porosity::porosity_n >;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      SortedArrayView< localIndex const > const & regionFilter,
                      FaceManager const & faceManager,
                      CompFlowAccessors const & compFlowAccessors,
                      MultiFluidAccessors const & multiFluidAccessors,
                      PorosityAccessors const & porosityAccessors,
                      real64 const dt,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_dt( dt ),
    m_regionFilter( regionFilter ),
    m_elemRegionList( faceManager.elementRegionList() ),
    m_elemSubRegionList( faceManager.elementSubRegionList() ),
    m_elemList( faceManager.elementList() ),
    m_volume( compFlowAccessors.get( fields::elementVolume {} ) ),
    m_porosity_n( porosityAccessors.get( fields::porosity::porosity_n {} ) ),
    m_totalDens_n( multiFluidAccessors.get( fields::multifluid::totalDensity_n {} ) )
  {}

  GEOS_HOST_DEVICE
  void computeMassNormalizer( localIndex const kf,
                              real64 & massNormalizer ) const
  {
    integer elemCounter = 0;

    for( integer k = 0; k < m_elemRegionList.size( 1 ); ++k )
    {
      localIndex const er  = m_elemRegionList[kf][k];
      localIndex const esr = m_elemSubRegionList[kf][k];
      localIndex const ei  = m_elemList[kf][k];
      bool const onBoundary = (er == -1 || esr == -1 || ei == -1);
      bool const isInTarget = m_regionFilter.contains( er );

      // if not on boundary, increment the normalizer
      if( !onBoundary && isInTarget )
      {
        massNormalizer += m_totalDens_n[er][esr][ei][0] * m_porosity_n[er][esr][ei][0] * m_volume[er][esr][ei];
        elemCounter++;
      }
    }
    massNormalizer /= elemCounter; // average mass in the adjacent cells at the previous converged time step
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const kf,
                            LinfStackVariables & stack ) const override
  {
    // if the face is adjacent to target region, compute the local values
    if( m_dofNumber[kf] >= 0 )
    {
      real64 massNormalizer = 0.0;
      computeMassNormalizer( kf, massNormalizer );

      // scaled residual to be in mass units (needed because element and face residuals are blended in a single norm)
      stack.localValue[0] += LvArray::math::abs( m_localResidual[stack.localRow] * m_dt ) / LvArray::math::max( m_minNormalizer, massNormalizer );
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const kf,
                          L2StackVariables & stack ) const override
  {
    // if the face is adjacent to target region, compute the local values
    if( m_dofNumber[kf] >= 0 )
    {
      real64 massNormalizer = 0;
      computeMassNormalizer( kf, massNormalizer );

      // scaled residual to be in mass units (needed because element and face residuals are blended in a single norm)
      real64 const valMass = m_localResidual[stack.localRow] * m_dt;
      stack.localValue[0] += valMass * valMass;
      stack.localNormalizer[0] += massNormalizer;
    }
  }


protected:

  /// Time step size
  real64 const m_dt;

  /// Filter to identify the target regions of the solver
  SortedArrayView< localIndex const > const m_regionFilter;

  /// Views on the maps face to elements
  arrayView2d< localIndex const > const m_elemRegionList;
  arrayView2d< localIndex const > const m_elemSubRegionList;
  arrayView2d< localIndex const > const m_elemList;

  /// View on the volume
  ElementViewConst< arrayView1d< real64 const > > const m_volume;

  /// View on porosity at the previous converged time step
  ElementViewConst< arrayView2d< real64 const > > const m_porosity_n;

  /// View on total mass/molar density at the previous converged time step
  ElementViewConst< arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > > const m_totalDens_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] regionFilter filter to identify the target regions of the solver
   * @param[in] solverName the name of the solver
   * @param[in] elemManager reference to the element region manager
   * @param[in] faceManager reference to the face manager
   * @param[in] dt time step size
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( physicsSolverBaseKernels::NormType const normType,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   SortedArrayView< localIndex const > const & regionFilter,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   FaceManager const & faceManager,
                   real64 const dt,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1],
                   real64 (& residualNormalizer)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = faceManager.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = faceManager.ghostRank();

    using kernelType = ResidualNormKernel;
    typename kernelType::CompFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::MultiFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PorosityAccessors poroAccessors( elemManager, solverName );

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank,
                               regionFilter, faceManager, flowAccessors, fluidAccessors, poroAccessors, dt, minNormalizer );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( faceManager.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( faceManager.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{

  template< typename POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & facePres,
          real64 const scalingFactor )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const iface )
    {
      if( ghostRank[iface] < 0 && dofNumber[iface] >= 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[iface] - rankOffset );
        {
          real64 const newFacePres = facePres[iface] + scalingFactor * localSolution[localRow];
          check.min( newFacePres >= 0.0 );
        }
      }
    } );
    return check.get();
  }

};

/******************************** PrecomputeKernel ********************************/

struct PrecomputeKernel
{

  template< typename IP_TYPE, integer NF >
  static void
  launch( localIndex const subRegionSize,
          localIndex const faceManagerSize,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView2d< real64 const > const & elemCenter,
          arrayView1d< real64 const > const & elemVolume,
          arrayView3d< real64 const > const & elemPerm,
          arrayView1d< real64 const > const & elemGravCoef,
          arrayView2d< localIndex const > const & elemToFaces,
          arrayView1d< real64 const > const & transMultiplier,
          real64 const & lengthTolerance,
          arrayView1d< RAJA::ReduceSum< serialReduce, real64 > > const & mimFaceGravCoefNumerator,
          arrayView1d< RAJA::ReduceSum< serialReduce, real64 > > const & mimFaceGravCoefDenominator,
          arrayView1d< real64 > const & mimFaceGravCoef )
  {
    forAll< serialPolicy >( subRegionSize, [=] ( localIndex const ei )
    {
      stackArray2d< real64, NF *NF > transMatrix( NF, NF );

      real64 const perm[ 3 ] = { elemPerm[ei][0][0], elemPerm[ei][0][1], elemPerm[ei][0][2] };

      IP_TYPE::template compute< NF >( nodePosition,
                                       transMultiplier,
                                       faceToNodes,
                                       elemToFaces[ei],
                                       elemCenter[ei],
                                       elemVolume[ei],
                                       perm,
                                       lengthTolerance,
                                       transMatrix );

      for( integer ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
      {
        mimFaceGravCoefNumerator[elemToFaces[ei][ifaceLoc]] += elemGravCoef[ei] * transMatrix[ifaceLoc][ifaceLoc];
        mimFaceGravCoefDenominator[elemToFaces[ei][ifaceLoc]] += transMatrix[ifaceLoc][ifaceLoc];
      }
    } );

    forAll< serialPolicy >( faceManagerSize, [=] ( localIndex const iface )
    {
      if( !isZero( mimFaceGravCoefDenominator[iface].get() ) )
      {
        mimFaceGravCoef[iface] = mimFaceGravCoefNumerator[iface].get() / mimFaceGravCoefDenominator[iface].get();
      }
    } );
  }
};

/******************************** Kernel switches ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void kernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernelLaunchSelectorFaceSwitch: type should be integral" );

  switch( value )
  {
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return;}
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return;}
    case 6:
    { lambda( std::integral_constant< T, 6 >() ); return;}
    default: GEOS_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

} // namespace internal

template< typename KERNELWRAPPER, typename INNER_PRODUCT, typename ... ARGS >
void simpleKernelLaunchSelector( localIndex numFacesInElem, ARGS && ... args )
{
  internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NUM_FACES )
  {
    KERNELWRAPPER::template launch< INNER_PRODUCT, NUM_FACES() >( std::forward< ARGS >( args )... );
  } );
}


template< typename KERNELWRAPPER, typename IP_TYPE, typename ... ARGS >
void kernelLaunchSelector( integer numFacesInElem, integer numComps, integer numPhases, ARGS && ... args )
{
  // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
  if( numPhases == 2 )
  {
    if( numComps == 2 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 2, 2, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 3 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 3, 2, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 4 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 4, 2, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 5 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 5, 2, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else
    {
      GEOS_ERROR( "Unsupported number of components: " << numComps );
    }
  }
  else if( numPhases == 3 )
  {
    if( numComps == 2 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 2, 3, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 3 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 3, 3, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 4 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 4, 3, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 5 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 5, 3, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else
    {
      GEOS_ERROR( "Unsupported number of components: " << numComps );
    }
  }
  else
  {
    GEOS_ERROR( "Unsupported number of phases: " << numPhases );
  }
}

} // namespace compositionalMultiphaseHybridFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
