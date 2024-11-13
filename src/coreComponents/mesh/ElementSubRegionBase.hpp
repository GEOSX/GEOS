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
 * @file ElementSubRegionBase.hpp
 */

#ifndef GEOS_MESH_ELEMENTSUBREGIONBASE_HPP_
#define GEOS_MESH_ELEMENTSUBREGIONBASE_HPP_

#include "mesh/ElementType.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

class NodeManager;
class FaceManager;
class MeshLevel;

namespace constitutive
{
class ConstitutiveBase;
}

/**
 * @class ElementSubRegionBase
 * Abstract class for a collection of mesh elements that
 * will be derived and specialized for cell elements,
 * face (fracture) elements, embedded surface (fracture)
 * elements, well elements, etc
 */
class ElementSubRegionBase : public ObjectManagerBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ElementSubRegionBase( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Destructor.
   */
  ~ElementSubRegionBase();

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  /**
   * @brief Calculate the geometric quantities for each element in the subregion.
   * @param[in] nodeManager the nodeManager (for geometrical info and connectivity involving nodes)
   * @param[in] faceManager the faceManager (for geometrical info and connectivity involving faces)
   *
   */
  virtual void calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & faceManager ) = 0;

  /**
   * @brief Link the connectivity maps of the subregion to the managers storing the mesh information.
   * @param[in] mesh the meshLevel object (single level only)
   *
   * In the derived classes, this function is used to passe a pointer to the nodeManager,
   * faceManager, and (if needed) edgeManager to, respectively, the node list, face list, and
   * edge list of the subregion.
   */
  virtual void setupRelatedObjectsInRelations( MeshLevel const & mesh ) = 0;

  /**
   * @brief Call ObjectManagerBase::fixUpDownMaps for the connectivity maps needed by
   *        the derived class (i.e., element-to-node map, element-to-face map, etc)
   * @param[in] clearIfUnmapped clearIfUnmapped
   */
  virtual void fixUpDownMaps( bool const clearIfUnmapped ) { GEOS_UNUSED_VAR( clearIfUnmapped ); }


  /**
   * @brief Set all "perElement" values for this subregion.
   * @param numNodesPerElement The number of nodes per elem in this subregion.
   * @param numEdgesPerElement The number of edges per elem in this subregion.
   * @param numFacesPerElement The number of faces per elem in this subregion.
   */
  virtual void resizePerElementValues( localIndex const numNodesPerElement,
                                       localIndex const numEdgesPerElement,
                                       localIndex const numFacesPerElement );

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get the number of nodes per element.
   * @return number of nodes per element
   */
  localIndex const & numNodesPerElement() const { return m_numNodesPerElement; }

  /**
   * @brief Get the number of nodes per element.
   * @param[in] k cell index
   * @return number of nodes per element
   */
  virtual localIndex numNodesPerElement( localIndex const k ) const { GEOS_UNUSED_VAR( k ); return m_numNodesPerElement; }

  /**
   * @brief Get the number of faces per element.
   * @return number of faces per element
   */
  localIndex const & numFacesPerElement() const { return m_numFacesPerElement; }

  /**
   * @brief Gets the number of edges per element.
   * @return input The number of edges per element in this ElementSubRegion.
   */
  localIndex const & numEdgesPerElement() const { return m_numEdgesPerElement; }

  /**
   * @brief Get the center of each element in this subregion.
   * @return an arrayView1d of const element centers
   */
  arrayView2d< real64 const > getElementCenter() const
  { return m_elementCenter; }

  /**
   * @copydoc getElementCenter() const
   */
  array2d< real64 > & getElementCenter()
  { return m_elementCenter; }

  /**
   * @brief Get the volume of each element in this subregion.
   * @return an arrayView1d of const element volumes
   */
  arrayView1d< real64 const > getElementVolume() const
  { return m_elementVolume; }

  /**
   * @brief Get the group in which the constitutive models of this subregion are registered.
   * @return a pointer to the const group in which the constitutive models are registered
   */
  dataRepository::Group const & getConstitutiveModels() const
  { return m_constitutiveModels; }

  /**
   * @copydoc getConstitutiveModels() const
   */
  dataRepository::Group & getConstitutiveModels()
  { return m_constitutiveModels; }

  /**
   * @brief Get a pointer to the constitutive model.
   * @tparam T The type of the constitutive model.
   * @param name The name of the constitutive model.
   * @return A pointer to the constitutive model.
   */
  template< typename T = constitutive::ConstitutiveBase >
  T const & getConstitutiveModel( string const & name ) const
  { return m_constitutiveModels.getGroup< T >( name ); }

  /**
   * @copydoc getConstitutiveModel( string const & ) const
   */
  template< typename T = constitutive::ConstitutiveBase >
  T & getConstitutiveModel( string const & name )
  { return m_constitutiveModels.getGroup< T >( name ); }

  /**
   * @brief Get the type of element in this subregion.
   * @return the type of element in this subregion
   */
  ElementType getElementType() const
  { return m_elementType; }

  /**
   * @brief Setter for m_elementType
   * @param elemType They type of element for this ElementSubRegion.
   */
  void setElementType( ElementType const elemType )
  { m_elementType = elemType; }


  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @return String key for the number of nodes per element in this subregion.
    static constexpr char const * numNodesPerElementString() { return "numNodesPerElement"; }
    /// @return String key for the element-to-node relation
    static constexpr char const * nodeListString() { return "nodeList"; }
    /// @return String key for the number of edges per element in this subregion.
    static constexpr char const * numEdgesPerElementString() { return "numEdgesPerElement"; }
    /// @return String key for the element-to-edge relation
    static constexpr char const * edgeListString() { return "edgeList"; }
    /// @return String key for the number of faces per element in this subregion.
    static constexpr char const * numFacesPerElementString() { return "numFacesPerElement"; }
    /// @return String key for the element-to-face relation
    static constexpr char const * faceListString() { return "faceList"; }
    /// @return String key for the member level field for the element center.
    static constexpr char const * elementCenterString() { return "elementCenter"; }
    /// @return String key for the member level field for the element volume.
    static constexpr char const * elementVolumeString() { return "elementVolume"; }
  };

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeyStruct
   */
  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    /// @return String key for the group in which the constitutive models of this subregion are registered.
    static constexpr auto constitutiveModelsString() { return "ConstitutiveModels"; }
  };

private:
  /// Group in which the constitutive models of this subregion are registered
  dataRepository::Group m_constitutiveModels;

protected:
  /// Number of nodes per element in this subregion.
  localIndex m_numNodesPerElement;

  /// Number of edges per element in this subregion.
  localIndex m_numEdgesPerElement;

  /// Number of faces per element in this subregion.
  localIndex m_numFacesPerElement;

  /// Member level field for the element center.
  array2d< real64 > m_elementCenter;

  /// Member level field for the element volume.
  array1d< real64 > m_elementVolume;

  /// Type of element in this subregion.
  ElementType m_elementType;

  /**
   * @brief Compute the center of each element in the subregion.
   * @tparam NODE_MAP Type of the element to node mapping.
   * @param toNodesRelation Element to node mapping
   * @param X Node positions.
   */
  template< class NODE_MAP >
  void calculateElementCenters( NODE_MAP const & toNodesRelation,
                                arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
  {
    arrayView2d< real64 > const & elementCenters = m_elementCenter;
    auto const e2n = toNodesRelation.toViewConst();

    forAll< parallelHostPolicy >( size(), [=]( localIndex const k )
    {
      LvArray::tensorOps::copy< 3 >( elementCenters[k], X[e2n( k, 0 )] );
      localIndex const numNodes = this->numNodesPerElement( k );
      for( localIndex a = 1; a < numNodes; ++a )
      {
        LvArray::tensorOps::add< 3 >( elementCenters[k], X[e2n( k, a )] );
      }

      LvArray::tensorOps::scale< 3 >( elementCenters[k], 1.0 / numNodes );
    } );
  }
};

} /* namespace geos */

#endif /* GEOS_MESH_ELEMENTSUBREGIONBASE_HPP_ */
