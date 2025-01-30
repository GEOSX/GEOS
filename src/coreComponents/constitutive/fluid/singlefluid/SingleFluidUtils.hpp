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
 * @file SingleFluidUtils.hpp
 */
#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SINGLEFLUIDUTILS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SINGLEFLUIDUTILS_HPP_

#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{

namespace constitutive
{


/**
 * @brief Helper struct used to represent a variable and its compositional derivatives
 * @tparam DIM number of dimensions
 */
template< typename T, int DIM, int USD, int USD_DC >
struct SingleFluidVarSlice
{
  GEOS_HOST_DEVICE
  SingleFluidVarSlice( internal::ArraySliceOrRef< T, DIM, USD > inputValue,
                       internal::ArraySliceOrRef< T, DIM+1, USD_DC > inputDerivs ):
    value( inputValue ),
    derivs( inputDerivs )
  {}

  internal::ArraySliceOrRef< T, DIM, USD > value;        /// variable value
  internal::ArraySliceOrRef< T, DIM + 1, USD_DC > derivs; /// derivative w.r.t. pressure, temperature, compositions

  using ValueType = internal::ArraySliceOrRef< T, DIM, USD >;
};

/**
 * @brief Struct holding views into fluid data, used to simplify parameter passing in kernel wrapper constructors.
 * @tparam NDIM number of dimensions
 * @tparam USD unit-stride-dim of primary property
 * @tparam USD_DC unit-stride-dim of derivatives
 */
template< typename T, int NDIM, int USD, int USD_DC >
struct SingleFluidVarView
{
  SingleFluidVarView() = default;

  GEOS_HOST_DEVICE
  SingleFluidVarView ( SingleFluidVarView const & src ):
    value( src.value ),
    derivs( src.derivs )
  {}

  GEOS_HOST_DEVICE
  SingleFluidVarView ( ArrayView< T, NDIM, USD > const & valueSrc,
                       ArrayView< T, NDIM + 1, USD_DC > const & derivsSrc ):
    value( valueSrc ),
    derivs( derivsSrc )
  {}

  ArrayView< T, NDIM, USD > value;        ///< View into property values
  ArrayView< T, NDIM + 1, USD_DC > derivs; ///< View into property derivatives w.r.t. pressure, temperature, compositions

  using SliceType = SingleFluidVarSlice< T, NDIM - 2, USD - 2, USD_DC - 2 >;

  using ValueType = ArrayView< T, NDIM, USD >;

  GEOS_HOST_DEVICE
  SliceType operator()( localIndex const k, localIndex const q ) const
  {
    return { value[k][q], derivs[k][q] };
  }
};

/**
 * @brief Struct holding views into fluid data, used to simplify parameter passing in kernel wrapper constructors.
 * @tparam NDIM number of dimensions
 * @tparam PERM unit-stride-dim of primary property
 * @tparam PERM_DC unit-stride-dim of derivatives
 */
template< typename T, int NDIM, typename PERM, typename PERM_DC >
struct SingleFluidVar
{
  Array< real64, NDIM, PERM > value;         ///< Property values
  Array< real64, NDIM + 1, PERM_DC > derivs; ///< Property derivatives w.r.t. pressure, temperature, compositions

  using ViewType = SingleFluidVarView< T, NDIM, getUSD< PERM >, getUSD< PERM_DC > >;
  using ViewTypeConst = SingleFluidVarView< T const, NDIM, getUSD< PERM >, getUSD< PERM_DC > >;

  using SliceType = typename ViewType::SliceType;
  using SliceTypeConst = typename ViewTypeConst::SliceType;

  using ValueType = Array< real64, NDIM, PERM >;
  template< int MAXSIZE >
  using StackValueType = StackArray< real64, NDIM, MAXSIZE, PERM >;
  using ViewValueType = typename ViewType::ValueType;
  using SliceValueType = typename SliceType::ValueType;

  ViewType toView()
  {
    return { value.toView(), derivs.toView() };
  }

  ViewTypeConst toViewConst() const
  {
    return { value.toViewConst(), derivs.toViewConst() };
  }
};


} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_SINGLEFLUIDUTILS_HPP_
