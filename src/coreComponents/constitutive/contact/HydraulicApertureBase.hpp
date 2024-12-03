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
 * @file HydraulicApertureBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTUREBASE_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTUREBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "functions/TableFunction.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"


namespace geos
{

namespace constitutive
{

/**
 * @class HydraulicApertureBase
 *
 * This class serves as the interface for implementing contact enforcement constitutive relations.
 * This does not include the actual enforcement algorithm, but only the constitutive relations that
 * govern the behavior of the contact. So things like penalty, or friction, or kinematic constraint.
 */
class HydraulicApertureBase : public ConstitutiveBase
{
public:

  /**
   * @brief The standard data repository constructor
   * @param name The name of the relation in the data repository
   * @param parent The name of the parent Group that holds this relation object.
   */
  HydraulicApertureBase( string const & name,
                         Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~HydraulicApertureBase() override;

protected:

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    /// string/key for aperture under zero normal stress
    static constexpr char const * apertureZeroString() { return "referenceAperture"; }
  };

  /// Reference hydraulic aperture. Aperture at zero normal stress
  real64 m_aperture0;  /// TODO: this will replace what is currently called defaultAperture.
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_HYDRAULICAPERTURETABLE_HPP_ */
