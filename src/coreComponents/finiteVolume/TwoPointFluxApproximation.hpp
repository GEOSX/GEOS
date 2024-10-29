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
 * @file TwoPointFluxApproximation.hpp
 */

#ifndef GEOS_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
#define GEOS_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_

#include "finiteVolume/FluxApproximationBase.hpp"

namespace geos
{

/**
 * @class TwoPointFluxApproximation
 *
 * Provides management of the interior stencil points when using a two-point flux approximation.
 */
class TwoPointFluxApproximation : public FluxApproximationBase
{
public:

  /**
   * @brief Static Factory Catalog Functions.
   * @return the catalog name
   */
  static string catalogName() { return "TwoPointFluxApproximation"; }

  TwoPointFluxApproximation() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the TwoPointFluxApproximation in the data repository
   * @param parent the parent group of this group.
   */
  TwoPointFluxApproximation( string const & name, dataRepository::Group * const parent );

  /**
   * @brief View keys.
   */
  struct viewKeyStruct : FluxApproximationBase::viewKeyStruct
  {
    /// @return The key for edfStencil
    static constexpr char const * edfmStencilString() { return "edfmStencil"; }
    /// @return The key for faceElementToCellStencil
    static constexpr char const * faceToCellStencilString() { return "faceElementToCellStencil"; }
    /// @return The key for the meanPermCoefficient
    static constexpr char const * meanPermCoefficientString() { return "meanPermCoefficient"; }
    /// @return The key for the usePEDFM flag
    static constexpr char const * usePEDFMString() { return "usePEDFM"; }
  };

private:
  virtual void registerCellStencil( Group & stencilGroup ) const override;

  virtual void computeFractureStencil( MeshLevel & mesh ) const override;

  virtual void computeCellStencil( MeshLevel & mesh ) const override;

  virtual void registerFractureStencil( Group & stencilGroup ) const override;

  virtual void addToFractureStencil( MeshLevel & mesh,
                                     string const & faceElementRegionName ) const override;

  virtual void registerBoundaryStencil( Group & stencilGroup,
                                        string const & setName ) const override;

  virtual void computeBoundaryStencil( MeshLevel & mesh,
                                       string const & setName,
                                       SortedArrayView< localIndex const > const & faceSet ) const override;

  virtual void registerAquiferStencil( Group & stencilGroup,
                                       string const & setName ) const override;

  virtual void computeAquiferStencil( DomainPartition & domain,
                                      MeshLevel & mesh ) const override;

  virtual void addEmbeddedFracturesToStencils( MeshLevel & mesh,
                                               string const & embeddedSurfaceRegionName ) const override;

  /**
   * @brief adds fracture-fracture connections to the surfaceElement stencil
   * @param mesh the mesh object
   * @param embeddedSurfaceRegionName name of the fracture region
   */
  void addFractureFractureConnectionsEDFM( MeshLevel & mesh,
                                           string const & embeddedSurfaceRegionName ) const;

  /**
   * @brief adds fracture-matrix connections to the edfm stencil
   * @param mesh the mesh object
   * @param embeddedSurfaceRegionName name of the fracture region
   */
  void addFractureMatrixConnectionsEDFM( MeshLevel & mesh,
                                         string const & embeddedSurfaceRegionName ) const;

  /**
   * @brief Add the connections within the DFM fracture elements themselves.
   * @param mesh The mesh object.
   * @param faceElementRegionName Name of the DFM fracture region.
   */
  void addFractureFractureConnectionsDFM( MeshLevel & mesh,
                                          string const & faceElementRegionName ) const;

  /**
   * @brief Add the connections between the DFM fracture elements and the matrix elements.
   * @param mesh The mesh object.
   * @param faceElementRegionName Name of the DFM fracture region.
   */
  void addFractureMatrixConnectionsDFM( MeshLevel & mesh,
                                        string const & faceElementRegionName ) const;

  /**
   * @brief Remove the expired connections between in the matrix elements
   * @param mesh The mesh object.
   * @param faceElementRegionName Name of the DFM fracture region.
   *
   * When a fracture propagates, some elements of the matrix do not connect each other:
   * the fracture is now in between. The connections between these matrix elements need to be cancelled.
   */
  void cleanMatrixMatrixConnectionsDFM( MeshLevel & mesh,
                                        string const & faceElementRegionName ) const;

  /**
   * @brief Initialize fields on the newly created elements of the fracture.
   * @param mesh The mesh object.
   * @param faceElementRegionName Name of the DFM fracture region.
   * @deprecated All of this initialization should be performed elsewhere.
   * It is just here because it was convenient, but it is not appropriate
   * to have physics based initialization in the flux approximator.
   */
  void initNewFractureFieldsDFM( MeshLevel & mesh,
                                 string const & faceElementRegionName ) const;

  /// mean permeability coefficient
  real64 m_meanPermCoefficient;
  /// flag to determine whether or not to use projection EDFM
  integer m_useProjectionEmbeddedFractureMethod;
};

}


#endif //GEOS_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
