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
 * @file FluxApproximationBase.hpp
 *
 */

#ifndef GEOS_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
#define GEOS_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_

#include "dataRepository/Group.hpp"
#include "CellElementStencilTPFA.hpp"
#include "SurfaceElementStencil.hpp"
#include "FaceElementToCellStencil.hpp"
#include "EmbeddedSurfaceToCellStencil.hpp"
#include "SurfaceElementStencil.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

/**
 * @brief Upwinding scheme.
 */
enum class UpwindingScheme : integer
{
  PPU,    ///< PPU upwinding
  C1PPU,  ///< C1-PPU upwinding from https://doi.org/10.1016/j.advwatres.2017.07.028
  IHU ///< IHU as in https://link.springer.com/content/pdf/10.1007/s10596-019-09835-6.pdf
};

/**
 * @brief Strings for upwinding scheme.
 */
ENUM_STRINGS( UpwindingScheme,
              "PPU",
              "C1PPU",
              "IHU" );

/**
 * @struct UpwindingParameters
 *
 * Structure to store upwinding parameters, such as upwinding scheme and related tolerances etc.
 */
struct UpwindingParameters
{
  /// PPU or C1-PPU or IHU
  UpwindingScheme upwindingScheme;

  /// C1-PPU smoothing tolerance
  //real64 epsC1PPU;
};

/**
 * @class FluxApproximationBase
 *
 * Base class for various flux approximation classes.
 * Stores the main and boundary stencils, construction is implemented in derived classes.
 * Main stencil is the one for cell-to-cell fluxes.
 * Boundary stencils are for Dirichlet boundary conditions
 */
class FluxApproximationBase : public dataRepository::Group
{
public:

  /// Alias for CatalogInterface, necessary declarations for factory instantiation of derived classes
  using CatalogInterface = dataRepository::CatalogInterface< FluxApproximationBase, string const &, Group * const >;
  /**
   * @brief Return the data type in the data repository.
   * @return the data type in the data repository
   */
  static typename CatalogInterface::CatalogType & getCatalog();

  FluxApproximationBase() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the FluxApproximationBase in the data repository
   * @param parent the parent group of this group.
   */
  FluxApproximationBase( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Extract stencil stored under the mesh group.
   * @tparam TYPE type of Stencil to get
   * @param mesh the mesh level object
   * @param name name of the stencil object
   * @return reference to the stencil
   */
  template< typename TYPE >
  TYPE const & getStencil( MeshLevel const & mesh, string const & name ) const;

  /**
   * @copydoc getStencil(MeshLevel const &, string const &) const
   */
  template< typename TYPE >
  TYPE & getStencil( MeshLevel & mesh, string const & name ) const;

  /**
   * @brief Call a user-provided function for each stencil.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] mesh the mesh level containing the stencils
   * @param[in] lambda The LAMBDA function
   */
  template< typename LAMBDA >
  void forAllStencils( MeshLevel const & mesh, LAMBDA && lambda ) const;

  /**
   * @brief Call a user-provided function for the each stencil according to the provided TYPE.
   * @tparam TYPE The type to be passed to forWrappers
   * @tparam TYPES Other types to be passed to forWrappers
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] mesh the mesh level containing the stencils
   * @param[in] lambda The LAMBDA function
   */
  template< typename TYPE, typename ... TYPES, typename LAMBDA >
  void forStencils( MeshLevel const & mesh, LAMBDA && lambda ) const;

  /**
   * @brief Add a new fracture stencil.
   * @param[in,out] mesh the mesh on which to add the fracture stencil
   * @param[in] faceElementRegionName the face element region name
   */
  virtual void addToFractureStencil( MeshLevel & mesh,
                                     string const & faceElementRegionName ) const = 0;

  /**
   * @brief Add a new embedded fracture stencil.
   * @param[in,out] mesh the mesh on which to add the fracture stencil
   * @param[in] embeddedSurfaceRegionName the embedded surface element region name
   */
  virtual void addEmbeddedFracturesToStencils( MeshLevel & mesh,
                                               string const & embeddedSurfaceRegionName ) const = 0;

  /**
   * @brief View keys.
   */
  struct viewKeyStruct
  {
    /// @return The key for fieldName
    static constexpr char const * fieldNameString() { return "fieldName"; }

    /// @return The key for coefficientName
    static constexpr char const * coeffNameString() { return "coefficientName"; }

    /// @return The key for targetRegions
    static constexpr char const * targetRegionsString() { return "targetRegions"; }

    /// @return The key for areaRelTol
    static constexpr char const * areaRelativeToleranceString() { return "areaRelTol"; }

    /// @return The key for transMultiplier
    static constexpr char const * transMultiplierString() { return "TransMultiplier"; }


    // Keys below are for wrappers registered on MeshLevel, not the current object

    /// @return The key for cellStencil
    static constexpr char const * cellStencilString() { return "cellStencil"; }

    /// @return The key for fractureStencil
    static constexpr char const * fractureStencilString() { return "fractureStencil"; }

    /// @return The key for upwindingScheme
    static constexpr char const * upwindingSchemeString() { return "upwindingScheme"; }

    /// @return The key for epsC1PPU
    //static constexpr char const * epsC1PPUString() { return "epsC1PPU"; }

  };

  /**
   * @brief Group keys.
   */
  struct groupKeyStruct
  {
    /// @return Key under which the top-level group for all FV stencils will be registered on MeshLevel
    static constexpr auto stencilMeshGroupString() { return "finiteVolumeStencils"; }
  };

  /**
   * @brief get the list of the target regions on a given mesh body.
   * @param[in] meshBodyName name of the meshBody
   * @return a list of the target regions on the meshBody
   */
  array1d< string > & targetRegions( string const & meshBodyName ) { return m_targetRegions[meshBodyName]; }

  /**
   * @brief set the name of the field.
   * @param name name of the field to be set.
   */
  void addFieldName( string const & name );

  /**
   * @brief set the name of the coefficient.
   * @param name name of the coefficient.
   */
  void setCoeffName( string const & name );

  /**
   * @brief get the upwinding parameters.
   * @return upwinding parameters structure.
   */
  UpwindingParameters const & upwindingParams() const { return m_upwindingParams; }

protected:

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  /**
   * @brief Register the wrapper for cell stencil on a mesh.
   * @param stencilGroup the group holding the stencil objects
   */
  virtual void registerCellStencil( Group & stencilGroup ) const = 0;

  /**
   * @brief Actual computation of the cell-to-cell stencil.
   * @param[inout] mesh the mesh on which to perform the computation.
   */
  virtual void computeCellStencil( MeshLevel & mesh ) const = 0;

  /**
   * @brief Actual Computation of the fracture related stencils.
   * @param[inout] mesh the mesh on which to perform the computation.
   *
   * Compute the stencils within the fracture itself,
   * but between the fracture and the matrix too.
   */
  virtual void computeFractureStencil( MeshLevel & mesh ) const = 0;

  /**
   * @brief Register the wrapper for fracture stencil on a mesh.
   * @param stencilGroup the group holding the stencil objects
   */
  virtual void registerFractureStencil( Group & stencilGroup ) const = 0;

  /**
   * @brief Register the wrapper for boundary face stencil on a mesh.
   * @param stencilGroup the group holding the stencil objects
   * @param setName the face set name (used as the wrapper name)
   */
  virtual void registerBoundaryStencil( Group & stencilGroup,
                                        string const & setName ) const = 0;

  /**
   * @brief Allocate and populate a stencil to be used in dirichlet boundary condition application
   * @param mesh the target mesh level
   * @param setName name of the face set, to be used as wrapper name for the produced stencil
   * @param faceSet set of face indices to use
   */
  virtual void computeBoundaryStencil( MeshLevel & mesh,
                                       string const & setName,
                                       SortedArrayView< localIndex const > const & faceSet ) const = 0;

  /**
   * @brief Register the wrapper for aquifer stencil on a mesh.
   * @param stencilGroup the group holding the stencil objects
   * @param setName the face set name (used as the wrapper name)
   */
  virtual void registerAquiferStencil( Group & stencilGroup,
                                       string const & setName ) const = 0;

  /**
   * @brief Allocate and populate a stencil to be used in aquifer boundary condition application
   * @param domain the domain partion
   * @param mesh the target mesh level
   */
  virtual void computeAquiferStencil( DomainPartition & domain,
                                      MeshLevel & mesh ) const = 0;


  /// name of the primary solution field
  array1d< string > m_fieldNames;

  /// name of the coefficient field
  string m_coeffName;

  /// names of target regions to build the stencil for
  map< string, array1d< string > > m_targetRegions;

  /// relative tolerance
  real64 m_areaRelTol;

  /// length scale of the mesh body
  real64 m_lengthScale;

  /// upwinding parameters
  UpwindingParameters m_upwindingParams;

};

template< typename TYPE >
TYPE const & FluxApproximationBase::getStencil( MeshLevel const & mesh, string const & name ) const
{
  Group const & stencilGroup = mesh.getGroup( groupKeyStruct::stencilMeshGroupString() ).getGroup( getName() );
  return stencilGroup.getReference< TYPE >( name );
}

template< typename TYPE >
TYPE & FluxApproximationBase::getStencil( MeshLevel & mesh, string const & name ) const
{
  Group & stencilGroup = mesh.getGroup( groupKeyStruct::stencilMeshGroupString() ).getGroup( getName() );
  return stencilGroup.getReference< TYPE >( name );
}

template< typename LAMBDA >
void FluxApproximationBase::forAllStencils( MeshLevel const & mesh, LAMBDA && lambda ) const
{
  //TODO remove dependence on CellElementStencilTPFA and SurfaceElementStencil
  forStencils< CellElementStencilTPFA,
               SurfaceElementStencil,
               EmbeddedSurfaceToCellStencil,
               FaceElementToCellStencil >( mesh, std::forward< LAMBDA >( lambda ) );
}

template< typename TYPE, typename ... TYPES, typename LAMBDA >
void FluxApproximationBase::forStencils( MeshLevel const & mesh, LAMBDA && lambda ) const
{
  Group const & stencilGroup = mesh.getGroup( groupKeyStruct::stencilMeshGroupString() ).getGroup( getName() );
  stencilGroup.forWrappers< TYPE, TYPES... >( [&] ( auto const & wrapper )
  {
    lambda( wrapper.reference() );
  } );
}

} // namespace geos

#endif //GEOS_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
