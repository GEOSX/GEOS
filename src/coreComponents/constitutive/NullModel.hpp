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

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVE_NULLMODEL_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVE_NULLMODEL_HPP_

#include "ConstitutiveBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @class NullModel A null constitutive relation
 */
class NullModel : public constitutive::ConstitutiveBase
{
public:

  /// @copydoc geos::dataRepository::Group::Group
  NullModel( string const & name,
             Group * const parent );

  /// Destrutor
  virtual ~NullModel();

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "NullModel";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * Empty struct to serve as a KernelWrapper for the constitutive model.
   */
  struct KernelWrapper
  {

    /**
     * @brief Get shear modulus
     * @param[in] k Element index.
     * @return the shear modulus of element k
     */
    GEOS_HOST_DEVICE
    virtual real64 getShearModulus( localIndex const k ) const
    {
      GEOS_UNUSED_VAR( k );
      GEOS_ERROR( "getShearModulus() not implemented for this model" );

      return 0;
    }
  };

  /**
   * @brief Create a kernel wrapper for this constitutive relation.
   * @param includeState Whether or not to include the state in the wrapper.
   * @return 0
   */
  KernelWrapper createKernelUpdates( bool const includeState = false ) const
  {
    GEOS_UNUSED_VAR( includeState );
    return KernelWrapper();
  }
};

} // constitutive
} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVE_NULLMODEL_HPP_ */
