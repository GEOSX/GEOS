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
 * @file InvariantDecompositions.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_INVARIANTDECOMPOSITIONS_HPP
#define GEOS_CONSTITUTIVE_SOLID_INVARIANTDECOMPOSITIONS_HPP

namespace geos
{

namespace constitutive
{

/**
 * @brief Namespace to collect two-invariant decomposition helper functions.
 *
 * Stress and strain invariants are typically used in constitutive modeling
 * to define a frame-indifferent response.
 *
 * @note
 * Care must be used in interpreting the decompositions, as Voight
 * notation uses a different convention for stress vs strain.  In particular,
 * the deviator (\f$\mathbf{n}\f$) is stored as a
 * stress-like quantity, such that it satisfies the relationship
 *        \f[
 *        \sigma = p \mathbf{I} + \sqrt{2/3} q \mathbf{n}
 *        \f]
 * when written in "unrolled" form.  It does not immediately satisfy
 *        \f[
 *        \epsilon = \frac{1}{3} \epsilon_v \mathbf{I} + \sqrt{3/2} \epsilon_d \mathbf{n}
 *        \f]
 * due to the engineering notation used for off-diagonal strains.  The stress and
 * strain recomposition functions account for this subtlety, so the user should
 * be safe when using these utilities directly.
 */
namespace twoInvariant
{

/**
 * @brief Perform two-invariant decomposition of strain tensors
 *
 * @param[in] strain Strain tensor in Voight notation
 * @param[out] volStrain Volumetric strain invariant
 * @param[out] devStrain Deviatoric strain invariant
 * @param[out] deviator Unit orientation (tensor) for deviatoric part in "stress" Voight notation
 */
GEOS_HOST_DEVICE
inline
void strainDecomposition( real64 const ( &strain )[6],
                          real64 & volStrain,
                          real64 & devStrain,
                          real64 ( & deviator )[6] )
{
  volStrain = strain[0] + strain[1] + strain[2];

  for( localIndex i=0; i<3; ++i )
  {
    deviator[i] = strain[i] - volStrain/3.;
    deviator[i+3] = strain[i+3]/2.; // divide by two ("stress voight")
  }

  devStrain = 0;
  for( localIndex i=0; i<3; ++i )
  {
    devStrain += deviator[i] * deviator[i];
    devStrain += 2 * deviator[i+3] * deviator[i+3];
  }
  devStrain = std::sqrt( devStrain );

  if( devStrain < 1e-12 )
  {
    for( localIndex i=0; i<6; ++i )
    {
      deviator[i] = 0;
    }
  }
  else
  {
    for( localIndex i=0; i<6; ++i )
    {
      deviator[i] /= devStrain;
    }
  }
  devStrain *= sqrt( 2./3. );
  return;
}


/**
 * @brief Perform two-invariant decomposition of stress tensors
 *
 * @param[in] stress Stress tensor in Voight notation
 * @param[out] volStress Volumetric stress invariant (= mean stress)
 * @param[out] devStress Deviatoric stress invariant (= von Mises stress)
 * @param[out] deviator Unit orientation (tensor) for deviatoric part in "stress" Voight notation
 */
GEOS_HOST_DEVICE
inline
void stressDecomposition( real64 const ( &stress )[6],
                          real64 & volStress,
                          real64 & devStress,
                          real64 ( & deviator )[6] )
{
  volStress = ( stress[0] + stress[1] + stress[2] ) / 3;

  for( localIndex i=0; i<3; ++i )
  {
    deviator[i] = stress[i] - volStress;
    deviator[i+3] = stress[i+3];
  }

  devStress = 0;
  for( localIndex i=0; i<3; ++i )
  {
    devStress += deviator[i] * deviator[i];
    devStress += 2 * deviator[i+3] * deviator[i+3];
  }
  devStress = std::sqrt( devStress );

  if( devStress < 1e-12 )
  {
    for( localIndex i=0; i<6; ++i )
    {
      deviator[i] = 0;
    }
  }
  else
  {
    for( localIndex i=0; i<6; ++i )
    {
      deviator[i] /= devStress;
    }
  }
  devStress *= sqrt( 3./2. );

  return;
}


/**
 * @brief Perform two-invariant recomposition of strain tensors
 *
 * @param[in] volStrain Volumetric strain invariant (scalar)
 * @param[in] devStrain Deviatoric strain invariant (scalar)
 * @param[in] deviator Unit orientation (tensor) for deviatoric part in "stress" Voight notation
 * @param[out] strain Strain tensor in Voight notation
 */
GEOS_HOST_DEVICE
inline
void strainRecomposition( real64 const & volStrain,
                          real64 const & devStrain,
                          real64 const ( &deviator )[6],
                          real64 ( & strain )[6] )
{
  real64 const tmp = sqrt( 1.5 )*devStrain;

  for( localIndex i=0; i<3; ++i )
  {
    strain[i]   = volStrain/3. + tmp * deviator[i];
    strain[i+3] = 2 * tmp * deviator[i+3]; // engineering strain
  }

  return;
}


/**
 * @brief Perform two-invariant recomposition of stress tensors
 *
 * @param[in] volStress Volumetric stress invariant (mean stress)
 * @param[in] devStress Deviatoric stress invariant (von Mises stress)
 * @param[in] deviator Unit orientation (tensor) for deviatoric part in "stress" Voight notation
 * @param[out] stress Stress tensor in Voight notation
 */
GEOS_HOST_DEVICE
inline
void stressRecomposition( real64 const & volStress,
                          real64 const & devStress,
                          real64 const ( &deviator )[6],
                          real64 ( & stress )[6] )
{
  real64 const tmp = sqrt( 2./3. )*devStress;

  for( localIndex i=0; i<3; ++i )
  {
    stress[i]   = volStress + tmp * deviator[i];
    stress[i+3] = tmp * deviator[i+3];
  }

  return;
}

} /* namespace twoInvariant */

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_INVARIANTDECOMPOSITIONS_HPP */
