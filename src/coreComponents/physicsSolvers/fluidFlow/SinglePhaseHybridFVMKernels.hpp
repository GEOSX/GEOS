/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_SINGLEPHASEHYBRIDFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_SINGLEPHASEHYBRIDFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "finiteVolume/InnerProduct.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace InnerProduct;

namespace SinglePhaseHybridFVMKernels
{


/******************************** FluxKernelHelper ********************************/

struct FluxKernelHelper
{

  /**
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] elemDens the density at this elenent's center
   * @param[in] dElemDens_dp the derivative of the density at this element's center
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[out] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[out] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   */
  static
  void ComputeOneSidedVolFluxes( arrayView1d< real64 const > const & facePres,
                                 arrayView1d< real64 const > const & dFacePres,
                                 arrayView1d< real64 const > const & faceGravDepth,
                                 arraySlice1d< localIndex const > const elemToFaces,
                                 real64 const & elemPres,
                                 real64 const & dElemPres,
                                 real64 const & elemGravDepth,
                                 real64 const & elemDens,
                                 real64 const & dElemDens_dp,
                                 stackArray2d< real64, MAX_NUM_FACES
                                               *MAX_NUM_FACES > const & transMatrix,
                                 stackArray1d< real64, MAX_NUM_FACES > & oneSidedVolFlux,
                                 stackArray1d< real64, MAX_NUM_FACES > & dOneSidedVolFlux_dp,
                                 stackArray2d< real64, MAX_NUM_FACES
                                               *MAX_NUM_FACES > & dOneSidedVolFlux_dfp );

  /**
   * @brief In a given element, collect the upwinded mobilities at this element's faces
   * @param[in] mesh the mesh object (single level only)
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] elemDofKey
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[inout] upwMobility the upwinded mobilities at this element's faces
   * @param[inout] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or
   * neighbor)
   * @param[inout] upwDofNumber  the dof number of the upwind pressure
   *
   * Note: because of the upwinding, this function requires non-local information
   */
  static
  void UpdateUpwindedCoefficients( MeshLevel const * const mesh,
                                   array2d< localIndex > const & elemRegionList,
                                   array2d< localIndex > const & elemSubRegionList,
                                   array2d< localIndex > const & elemList,
                                   SortedArray< localIndex > const & regionFilter,
                                   arraySlice1d< localIndex const > const elemToFaces,
                                   ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & mob,
                                   ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & dMob_dp,
                                   localIndex const er,
                                   localIndex const esr,
                                   localIndex const ei,
                                   globalIndex const elemDofNumber,
                                   string const elemDofKey,
                                   stackArray1d< real64, MAX_NUM_FACES > const & oneSidedVolFlux,
                                   stackArray1d< real64, MAX_NUM_FACES > & upwMobility,
                                   stackArray1d< real64, MAX_NUM_FACES > & dUpwMobility_dp,
                                   stackArray1d< globalIndex, MAX_NUM_FACES > & upwDofNumber );

  /**
   * @brief In a given element, assemble the mass conservation equation
   * @param[in] dt the time step size
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] upwMobility the upwinded mobilities at this element's faces
   * @param[in] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or
   * neighbor)
   * @param[in] upwDofNumber  the dof number of the upwind pressure
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  static
  void AssembleOneSidedMassFluxes( real64 const & dt,
                                   arrayView1d< globalIndex const > const & faceDofNumber,
                                   arraySlice1d< localIndex const > const elemToFaces,
                                   globalIndex const elemDofNumber,
                                   stackArray1d< real64, MAX_NUM_FACES > const & oneSidedVolFlux,
                                   stackArray1d< real64, MAX_NUM_FACES > const & dOneSidedVolFlux_dp,
                                   stackArray2d< real64, MAX_NUM_FACES
                                                 *MAX_NUM_FACES > const & dOneSidedVolFlux_dfp,
                                   stackArray1d< real64, MAX_NUM_FACES > const & upwMobility,
                                   stackArray1d< real64, MAX_NUM_FACES > const & dUpwMobility_dp,
                                   stackArray1d< globalIndex, MAX_NUM_FACES > const & upwDofNumber,
                                   ParallelMatrix * const matrix,
                                   ParallelVector * const rhs );


  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  static
  void AssembleConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                            arraySlice1d< localIndex const > const elemToFaces,
                            globalIndex const elemDofNumber,
                            stackArray1d< real64, MAX_NUM_FACES > const & oneSidedVolFlux,
                            stackArray1d< real64, MAX_NUM_FACES > const & dOneSidedVolFlux_dp,
                            stackArray2d< real64, MAX_NUM_FACES
                                          *MAX_NUM_FACES > const & dOneSidedVolFlux_dfp,
                            ParallelMatrix * const matrix,
                            ParallelVector * const rhs );



};

/******************************** FluxKernel ********************************/

template< typename REGIONTYPE >
struct FluxKernel
{};

template<>
struct FluxKernel< CellElementSubRegion >
{

  inline static void
  Compute( localIndex const er,
           localIndex const esr,
           localIndex const ei,
           SortedArray< localIndex > const & regionFilter,
           MeshLevel const * const mesh,
           array2d< localIndex > const & elemRegionList,
           array2d< localIndex > const & elemSubRegionList,
           array2d< localIndex > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & dFacePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arraySlice1d< localIndex const > const elemToFaces,
           globalIndex const elemDofNumber,
           string const elemDofKey,
           real64 const & elemPres,
           real64 const & dElemPres,
           real64 const & elemGravCoef,
           real64 const & elemDens,
           real64 const & dElemDens_dp,
           ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & mobility,
           ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & dMobility_dPres,
           stackArray2d< real64, MAX_NUM_FACES
                         *MAX_NUM_FACES > const & transMatrix,
           real64 const & dt,
           ParallelMatrix * const matrix,
           ParallelVector * const rhs )
  {
    // max number of faces allowed in an element
    localIndex constexpr maxNumFaces = MAX_NUM_FACES;

    localIndex const numFacesInElem = elemToFaces.size();

    // one sided flux
    stackArray1d< real64, maxNumFaces > oneSidedVolFlux( numFacesInElem );
    stackArray1d< real64, maxNumFaces > dOneSidedVolFlux_dp( numFacesInElem );
    stackArray2d< real64, maxNumFaces *maxNumFaces > dOneSidedVolFlux_dfp( numFacesInElem, numFacesInElem );

    // upwinded mobility
    stackArray1d< real64, maxNumFaces > upwMobility( numFacesInElem );
    stackArray1d< real64, maxNumFaces > dUpwMobility_dp( numFacesInElem );
    stackArray1d< globalIndex, maxNumFaces > upwDofNumber( numFacesInElem );

    /*
     * compute auxiliary quantities at the one sided faces of this element:
     * 1) One-sided volumetric fluxes
     * 2) Upwinded mobilities
     */

    // for each one-sided face of the elem,
    // compute the volumetric flux using transMatrix
    FluxKernelHelper::ComputeOneSidedVolFluxes( facePres,
                                                dFacePres,
                                                faceGravCoef,
                                                elemToFaces,
                                                elemPres,
                                                dElemPres,
                                                elemGravCoef,
                                                elemDens,
                                                dElemDens_dp,
                                                transMatrix,
                                                oneSidedVolFlux,
                                                dOneSidedVolFlux_dp,
                                                dOneSidedVolFlux_dfp );

    // at this point, we know the local flow direction in the element
    // so we can upwind the transport coefficients (mobilities) at the one sided faces
    // ** this function needs non-local information **
    FluxKernelHelper::UpdateUpwindedCoefficients( mesh,
                                                  elemRegionList,
                                                  elemSubRegionList,
                                                  elemList,
                                                  regionFilter,
                                                  elemToFaces,
                                                  mobility,
                                                  dMobility_dPres,
                                                  er, esr, ei,
                                                  elemDofNumber,
                                                  elemDofKey,
                                                  oneSidedVolFlux,
                                                  upwMobility,
                                                  dUpwMobility_dp,
                                                  upwDofNumber );

    /*
     * perform assembly in this element in two steps:
     * 1) mass conservation equations
     * 2) face constraints
     */

    // use the computed one sided vol fluxes and the upwinded mobilities
    // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
    FluxKernelHelper::AssembleOneSidedMassFluxes( dt,
                                                  faceDofNumber,
                                                  elemToFaces,
                                                  elemDofNumber,
                                                  oneSidedVolFlux,
                                                  dOneSidedVolFlux_dp,
                                                  dOneSidedVolFlux_dfp,
                                                  upwMobility,
                                                  dUpwMobility_dp,
                                                  upwDofNumber,
                                                  matrix,
                                                  rhs );

    // use the computed one sided vol fluxes to assemble the constraints
    // enforcing flux continuity at this element's faces
    FluxKernelHelper::AssembleConstraints( faceDofNumber,
                                           elemToFaces,
                                           elemDofNumber,
                                           oneSidedVolFlux,
                                           dOneSidedVolFlux_dp,
                                           dOneSidedVolFlux_dfp,
                                           matrix,
                                           rhs );

  }

};



} // namespace SinglePhaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_SINGLEPHASEHYBRIDFVMKERNELS_HPP
