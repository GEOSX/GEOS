/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BartonBandis.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_BARTONBANDIS_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_BARTONBANDIS_HPP_

#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class BartonBandisUpdates
 *
 * This class is used for in-kernel contact relation updates
 */
class BartonBandisUpdates
{
public:

  BartonBandisUpdates( real64 const aperture0,
                       real64 const referenceNormalStress )
    : m_aperture0( aperture0 ),
    m_referenceNormalStress( referenceNormalStress )
  {}

  /// Default copy constructor
  BartonBandisUpdates( BartonBandisUpdates const & ) = default;

  /// Default move constructor
  BartonBandisUpdates( BartonBandisUpdates && ) = default;

  /// Deleted default constructor
  BartonBandisUpdates() = default;

  /// Deleted copy assignment operator
  BartonBandisUpdates & operator=( BartonBandisUpdates const & ) = delete;

  /// Deleted move assignment operator
  BartonBandisUpdates & operator=( BartonBandisUpdates && ) =  delete;

  /**
   * @brief Evaluate the effective aperture, and its derivative wrt aperture
   * @param[in] aperture the model aperture/gap
   * @param[out] dHydraulicAperture_dAperture the derivative of the effective aperture wrt aperture
   * @return The hydraulic aperture that is always > 0
   */
  GEOS_HOST_DEVICE
  real64 computeHydraulicAperture( real64 const aperture,
                                   real64 const normalTraction,
                                   real64 & dHydraulicAperture_aperture,
                                   real64 & dHydraulicAperture_dNormalStress ) const;

private:
  real64 m_aperture0;

  real64 m_referenceNormalStress;
};


/**
 * @class BartonBandis
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class BartonBandis : public HydraulicApertureBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  BartonBandis( string const & name,
                Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~BartonBandis() override;


  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BartonBandisUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  struct viewKeyStruct : public HydraulicApertureBase::viewKeyStruct
  {
    /// string/key for reference normal stress
    static constexpr char const * referenceNormalStressString() { return "referenceNormalStress"; }
  };
  /// Reference normal stress
  real64 m_referenceNormalStress;
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 BartonBandisUpdates::computeHydraulicAperture( real64 const aperture,
                                                      real64 const normalTraction,
                                                      real64 & dHydraulicAperture_aperture,
                                                      real64 & dHydraulicAperture_dNormalStress ) const
{
  real64 const hydraulicAperture = ( aperture >= 0.0 ) ? (aperture + m_aperture0) : m_aperture0 / ( 1 + 9*normalTraction/m_referenceNormalStress );
  dHydraulicAperture_dNormalStress = ( aperture >= 0.0 ) ? 0.0 : -hydraulicAperture / ( 1 + 9*normalTraction/m_referenceNormalStress ) * 9/m_referenceNormalStress;
  dHydraulicAperture_aperture = ( aperture >= 0.0 ) ? 1.0 : 0.0;

  return hydraulicAperture; ///It would be nice to change this to return a tuple.
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTURETABLE_HPP_ */
