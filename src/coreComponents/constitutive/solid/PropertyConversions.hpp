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
 * @file PropertyConversions.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP
#define GEOS_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP

namespace geos
{

namespace constitutive
{

/// @namespace Namespace to collect common property conversion functions (elastic, poroelastic, etc.)
namespace conversions
{

/// @namespace Bulk modulus and shear modulus as input
namespace bulkModAndShearMod
{

/**
 * @brief Compute Young's modulus
 * @param[in] K Bulk modulus
 * @param[in] G Shear modulus
 * @return Young's modulus
 */
GEOS_HOST_DEVICE
inline
real64 toYoungMod( real64 const & K, real64 const & G )
{
  return 9 * K * G / ( 3 * K + G );
}

/**
 * @brief Compute Poisson's ratio
 * @param[in] K Bulk modulus
 * @param[in] G Shear modulus
 * @return Poisson's ratio
 */
GEOS_HOST_DEVICE
inline
real64 toPoissonRatio( real64 const & K, real64 const & G )
{
  return ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );
}

/**
 * @brief Compute First Lamé parameter
 * @param[in] K Bulk modulus
 * @param[in] G Shear modulus
 * @return First Lamé parameter
 */
GEOS_HOST_DEVICE
inline
real64 toFirstLame( real64 const & K, real64 const & G )
{
  return K - 2 * G / 3;
}

} /* namespace bulkModeAndShearMod */

/// @namespace Young's modulus and Poisson ratio as input
namespace youngModAndPoissonRatio
{

/**
 * @brief Compute bulk modulus
 * @param[in] E Young's modulus
 * @param[in] nu Poisson's ratio
 * @return Bulk modulus
 */
GEOS_HOST_DEVICE
inline
real64 toBulkMod( real64 const & E, real64 const & nu )
{
  return E / (3 * ( 1 - 2*nu ) );
}

/**
 * @brief Compute bulk modulus
 * @param[in] E Young's modulus
 * @param[in] nu Poisson's ratio
 * @return Shear modulus
 */
GEOS_HOST_DEVICE
inline
real64 toShearMod( real64 const & E, real64 const & nu )
{
  return E / (2 * ( 1 + nu ) );
}

} /* namespace youngModAndPoissonRatio*/

/// @namespace Shear modulus and Poisson's ratio as input
namespace shearModAndPoissonRatio
{

/**
 * @brief Compute bulk modulus
 * @param[in] G Shear modulus
 * @param[in] nu Poisson's ratio
 * @return Bulk modulus
 */
GEOS_HOST_DEVICE
inline
real64 toBulkMod( real64 const & G, real64 const & nu )
{
  return 2 * G * ( 1 + nu ) / ( 3 * ( 1 - 2 * nu ) );
}

/**
 * @brief Compute Young's modulus
 * @param[in] G Shear modulus
 * @param[in] nu Poisson's ratio
 * @return Young's modulus
 */
GEOS_HOST_DEVICE
inline
real64 toYoungMod( real64 const & G, real64 const & nu )
{
  return 2 * G * ( 1 + nu );
}

} /* namespace shearModAndPoissonRatio*/

/// @namespace Bulk modulus and Poisson's ratio as input
namespace bulkModAndPoissonRatio
{

/**
 * @brief Compute Young's modulus
 * @param[in] K Bulk modulus
 * @param[in] nu Poisson's ratio
 * @return Young's modulus
 */
GEOS_HOST_DEVICE
inline
real64 toYoungMod( real64 const & K, real64 const & nu )
{
  return 3 * K * ( 1 - 2 * nu );
}

/**
 * @brief Compute Shear modulus
 * @param[in] K Bulk modulus
 * @param[in] nu Poisson's ratio
 * @return Shear modulus
 */
GEOS_HOST_DEVICE
inline
real64 toShearMod( real64 const & K, real64 const & nu )
{
  return 3 * K * ( 1 - 2 * nu) / ( 2 * ( 1 + nu ) );
}

} /* namespace bulkModAndPoissonRatio */

/// @namespace Bulk modulus and Young's modulus
namespace bulkModAndYoungMod
{

/**
 * @brief Compute Shear modulus
 * @param[in] K Bulk modulus
 * @param[in] E Young's ratio
 * @return Shear modulus
 */
GEOS_HOST_DEVICE
inline
real64 toShearMod( real64 const & K, real64 const & E )
{
  return 3 * K * E / ( 9 * K - E );
}

/**
 * @brief Compute Poisson ratio
 * @param[in] K Bulk modulus
 * @param[in] E Young's modulus
 * @return Poisson's ratio
 */
GEOS_HOST_DEVICE
inline
real64 toPoissonRatio( real64 const & K, real64 const & E )
{
  return ( 3 * K - E ) / ( 6 * K);
}

} /* namespace bulkModAndYoungMod */

/// @namespace Shear modulus and Young's modulus
namespace shearModAndYoungMod
{
/**
 * @brief Compute Poisson ratio
 * @param[in] G Shear modulus
 * @param[in] E Young's modulus
 * @return Poisson's ratio
 */
GEOS_HOST_DEVICE
inline
real64 toPoissonRatio( real64 const & G, real64 const & E )
{
  return 0.5 * E / G - 1.0;
}

/**
 * @brief Compute Bulk modulus
 * @param[in] G Shear modulus
 * @param[in] E Young's modulus
 * @return Bulk modulus
 */
GEOS_HOST_DEVICE
inline
real64 toBulkMod( real64 const & G, real64 const & E )
{
  return E * G / ( 3 * ( 3 * G - E ) );
}

} /* namespace shearModAndYoungMod*/

} /* namespace conversions */

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP */
