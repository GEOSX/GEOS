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
 * @file PressurePorosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_POROSITY_PROPPANTPOROSITY_HPP_
#define GEOS_CONSTITUTIVE_POROSITY_PROPPANTPOROSITY_HPP_

#include "PorosityBase.hpp"

namespace geos
{
namespace constitutive
{

class ProppantPorosityUpdates : public PorosityBaseUpdates
{
public:

  ProppantPorosityUpdates( arrayView2d< real64 > const & newPorosity,
                           arrayView2d< real64 const > const & porosity_n,
                           arrayView2d< real64 > const & dPorosity_dPressure,
                           arrayView2d< real64 > const & dPorosity_dTemperature,
                           arrayView2d< real64 const > const & initialPorosity,
                           arrayView1d< real64 const > const & referencePorosity,
                           real64 const & maxProppantConcentration ):
    PorosityBaseUpdates( newPorosity,
                         porosity_n,
                         dPorosity_dPressure,
                         dPorosity_dTemperature,
                         initialPorosity,
                         referencePorosity ),
    m_maxProppantConcentration( maxProppantConcentration )
  {}

  GEOS_HOST_DEVICE
  void computePorosity( real64 const & proppantPackVolumeFraction,
                        real64 & porosity ) const
  {
    porosity = 1.0 - m_maxProppantConcentration * proppantPackVolumeFraction;
  }

  GEOS_HOST_DEVICE
  void updateFromProppantVolumeFraction( localIndex const k,
                                         localIndex const q,
                                         real64 const & proppantPackVolumeFraction ) const
  {
    computePorosity( proppantPackVolumeFraction,
                     m_newPorosity[k][q] );
  }

private:

  real64 const m_maxProppantConcentration;

};


class ProppantPorosity : public PorosityBase
{
public:
  ProppantPorosity( string const & name, Group * const parent );

  virtual ~ProppantPorosity() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ProppantPorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public PorosityBase::viewKeyStruct
  {
    static constexpr char const * maxProppantConcentrationString() { return "maxProppantConcentration"; }
  } viewKeys;

  using KernelWrapper = ProppantPorosityUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates() const
  {
    return KernelWrapper( m_newPorosity,
                          m_porosity_n,
                          m_dPorosity_dPressure,
                          m_dPorosity_dTemperature,
                          m_initialPorosity,
                          m_referencePorosity,
                          m_maxProppantConcentration );
  }


private:
  virtual void postInputInitialization() override;

  real64 m_maxProppantConcentration;


};


}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_POROSITY_PROPPANTPOROSITY_HPP_
