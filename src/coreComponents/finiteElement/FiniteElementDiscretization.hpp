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
 * @file FiniteElementSpace.hpp
 */

#ifndef GEOS_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_
#define GEOS_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "mesh/ElementType.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "FiniteElementDispatch.hpp"


namespace geos
{

// TODO remove when these quantities are placed inside the FiniteElementBase
// class.
namespace dataRepository
{
namespace keys
{
string const dNdX = "dNdX";
string const detJ = "detJ";
}
}



class FiniteElementDiscretization : public dataRepository::Group
{
public:

  /// Enumerator of available interpolation types
  enum class Formulation : integer
  {
    Default,
    SEM,
    DG
  };


  FiniteElementDiscretization() = delete;

  explicit FiniteElementDiscretization( string const & name, Group * const parent );

  ~FiniteElementDiscretization() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  static string catalogName() { return "FiniteElementSpace"; }

  ///@}

  template< typename SUBREGION_TYPE,
            typename FE_TYPE >
  void calculateShapeFunctionGradients( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                        SUBREGION_TYPE * const elementSubRegion,
                                        typename FE_TYPE::template MeshData< SUBREGION_TYPE > meshData,
                                        FE_TYPE & fe ) const;


  /**
   * @brief Factory method to instantiate a type of finite element formulation.
   * @param parentElementShape String key that indicates the type of
   *   element/basis/formulation that should be instantiated.
   * @return A unique_ptr< FinteElementBase > which contains the new
   *   instantiation.
   */
  std::unique_ptr< finiteElement::FiniteElementBase >
  factory( ElementType const parentElementShape ) const;

  int getOrder() const { return m_order; }

private:

  struct viewKeyStruct
  {
    static constexpr char const * orderString() { return "order"; }
    static constexpr char const * formulationString() { return "formulation"; }
    static constexpr char const * useVemString() { return "useVirtualElements"; }
  };

  /// The order of the finite element basis
  int m_order;

  /// Optional string indicating any specialized formulation type.
  Formulation m_formulation;

  /// Optional parameter indicating if the class should use Virtual Elements.
  int m_useVem;

  void postInputInitialization() override final;

};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( FiniteElementDiscretization::Formulation,
              "default",
              "SEM",
              "DG" );

template< typename SUBREGION_TYPE,
          typename FE_TYPE >
void
FiniteElementDiscretization::
  calculateShapeFunctionGradients( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                   SUBREGION_TYPE * const elementSubRegion,
                                   typename FE_TYPE::template MeshData< SUBREGION_TYPE > meshData,
                                   FE_TYPE & finiteElement ) const
{
  GEOS_MARK_FUNCTION;

  // do not precompute shape functions in case of SEM or DG formulation (not needed)
  if( m_formulation == Formulation::SEM || m_formulation == Formulation::DG )
    return;

  array4d< real64 > & dNdX = elementSubRegion->dNdX();
  array2d< real64 > & detJ = elementSubRegion->detJ();
  auto const & elemsToNodes = elementSubRegion->nodeList().toViewConst();

  constexpr localIndex numNodesPerElem = FE_TYPE::maxSupportPoints;
  constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
  dNdX.resizeWithoutInitializationOrDestruction( elementSubRegion->size(), numQuadraturePointsPerElem, numNodesPerElem, 3 );
  detJ.resize( elementSubRegion->size(), numQuadraturePointsPerElem );

  finiteElement.setGradNView( dNdX.toViewConst() );
  finiteElement.setDetJView( detJ.toViewConst() );

  for( localIndex k = 0; k < elementSubRegion->size(); ++k )
  {
    typename FE_TYPE::StackVariables feStack;
    finiteElement.template setup< FE_TYPE >( k, meshData, feStack );
    real64 xLocal[numNodesPerElem][3];
    localIndex numSupportPoints = finiteElement.template numSupportPoints< FE_TYPE >( feStack );
    for( localIndex a=0; a< numSupportPoints; ++a )
    {
      localIndex const nodeIndex = elemsToNodes[ k][ a ];
      for( int i=0; i<3; ++i )
      {
        xLocal[ a ][ i ] = X[ nodeIndex ][ i ];
      }
    }



    for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
    {
      real64 dNdXLocal[numNodesPerElem][3];
      detJ( k, q ) = finiteElement.calcGradN( q, xLocal, feStack, dNdXLocal );

      for( localIndex b = 0; b < numSupportPoints; ++b )
      {
        LvArray::tensorOps::copy< 3 >( dNdX[ k ][ q ][ b ], dNdXLocal[b] );
      }
    }
  }

}


} /* namespace geos */

#endif /* GEOS_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_ */
