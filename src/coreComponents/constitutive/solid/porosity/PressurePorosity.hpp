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
 * @file PressurePorosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_POROSITY_PRESSUREPOROSITY_HPP_
#define GEOS_CONSTITUTIVE_POROSITY_PRESSUREPOROSITY_HPP_

#include "PorosityBase.hpp"

namespace geos
{
namespace constitutive
{

class PressurePorosityUpdates : public PorosityBaseUpdates
{
public:

  PressurePorosityUpdates( arrayView2d< real64 > const & newPorosity,
                           arrayView2d< real64 const > const & porosity_n,
                           arrayView2d< real64 > const & dPorosity_dPressure,
                           arrayView2d< real64 > const & dPorosity_dTemperature,
                           arrayView2d< real64 const > const & initialPorosity,
                           arrayView1d< real64 const > const & referencePorosity,
                           real64 const & referencePressure,
                           real64 const & compressibility ):
    PorosityBaseUpdates( newPorosity,
                         porosity_n,
                         dPorosity_dPressure,
                         dPorosity_dTemperature,
                         initialPorosity,
                         referencePorosity ),
    m_referencePressure( referencePressure ),
    m_compressibility( compressibility )
  {}

  GEOS_HOST_DEVICE
  void computePorosity( real64 const & pressure,
                        real64 const & temperature,
                        real64 & porosity,
                        real64 & dPorosity_dPressure,
                        real64 & dPorosity_dTemperature,
                        real64 const & referencePorosity ) const
  {
    GEOS_UNUSED_VAR( temperature );

    // TODO use full exponential.
//    porosity            =  referencePorosity * exp( m_compressibility * (pressure - m_referencePressure) );
//    dPorosity_dPressure =  m_compressibility * porosity;
    porosity = referencePorosity * ( m_compressibility * (pressure - m_referencePressure) + 1 );
    dPorosity_dPressure = m_compressibility * referencePorosity;
    dPorosity_dTemperature = 0.0;
  }

  GEOS_HOST_DEVICE
  void updateFromPressureAndTemperature( localIndex const k,
                                         localIndex const q,
                                         real64 const & pressure,
                                         real64 const & temperature ) const
  {
    computePorosity( pressure,
                     temperature,
                     m_newPorosity[k][q],
                     m_dPorosity_dPressure[k][q],
                     m_dPorosity_dTemperature[k][q],
                     m_referencePorosity[k] );
  }

private:

  /// Reference pressure used in the porosity model
  real64 const m_referencePressure;

  /// Compressibility used in the porosity model
  real64 const m_compressibility;
};


class PressurePorosity : public PorosityBase
{
public:
  PressurePorosity( string const & name, Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PressurePorosity"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public PorosityBase::viewKeyStruct
  {
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * compressibilityString() { return "compressibility"; }
  } viewKeys;

  using KernelWrapper = PressurePorosityUpdates;

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
                          m_referencePressure,
                          m_compressibility );
  }


private:
  virtual void postInputInitialization() override;

  real64 m_referencePressure;

  real64 m_compressibility;

};


}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_POROSITY_PRESSUREPOROSITY_HPP_
