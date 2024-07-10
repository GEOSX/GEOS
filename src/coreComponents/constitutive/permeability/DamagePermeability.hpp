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
 * @file DamagePermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_DAMAGEPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_DAMAGEPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

class DamagePermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  DamagePermeabilityUpdate( arrayView3d< real64 > const & permeability,
                            arrayView3d< real64 > const & dPerm_dPressure,
                            real64 const & bulkPermeability )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_bulkPermeability( bulkPermeability )
  {}

  GEOS_HOST_DEVICE
  void updateDamagePermeability ( localIndex const k,
                                  real64 const & damage ) const
  {
    real64 const matrixPermeability = m_bulkPermeability*exp( 7.0*damage );

    for( localIndex dim=0; dim<3; ++dim )
    {
      m_permeability[k][0][dim] = matrixPermeability;
    }
  }

private:

  real64 m_bulkPermeability;

};


class DamagePermeability : public PermeabilityBase
{
public:

  DamagePermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "DamagePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = DamagePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_bulkPermeability );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * bulkPermeabilityString() { return "bulkPermeability"; }
  };

protected:

  virtual void postInputInitialization() override;

private:

  real64 m_bulkPermeability;

};

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_DAMAGEPERMEABILITY_HPP_
