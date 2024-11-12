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
 * @file BrooksCoreyRelativePermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
#define GEOS_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"

namespace geos
{
namespace constitutive
{

class BrooksCoreyRelativePermeabilityUpdate final : public RelativePermeabilityBaseUpdate
{
public:

  BrooksCoreyRelativePermeabilityUpdate( arrayView2d< real64 const > const & phaseMinVolumeFraction,
                                         arrayView2d< real64 const > const & phaseRelPermExponent,
                                         arrayView2d< real64 const > const & phaseRelPermMaxValue,
                                         arrayView1d< real64 const > const & volFracScale,
                                         arrayView1d< integer const > const & phaseTypes,
                                         arrayView1d< integer const > const & phaseOrder,
                                         arrayView4d< real64, constitutive::relperm::USD_RELPERM > const & phaseRelPerm,
                                         arrayView5d< real64, constitutive::relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
                                         arrayView3d< real64, constitutive::relperm::USD_PHASE > const & phaseTrappedVolFrac )
    : RelativePermeabilityBaseUpdate( phaseTypes,
                                      phaseOrder,
                                      phaseRelPerm,
                                      dPhaseRelPerm_dPhaseVolFrac,
                                      phaseTrappedVolFrac ),
    m_phaseMinVolumeFraction( phaseMinVolumeFraction ),
    m_phaseRelPermExponent( phaseRelPermExponent ),
    m_phaseRelPermMaxValue( phaseRelPermMaxValue ),
    m_volFracScale( volFracScale )
  {}

  GEOS_HOST_DEVICE
  void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                arraySlice1d< real64, constitutive::relperm::USD_PHASE - 2 > const & phaseTrappedVolFrac,
                arraySlice2d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                arraySlice3d< real64, constitutive::relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override
  {
    compute( phaseVolFraction,
             m_phaseTrappedVolFrac[k][q],
             m_phaseRelPerm[k][q],
             m_dPhaseRelPerm_dPhaseVolFrac[k][q] );
  }

private:

  arrayView2d< real64 const > m_phaseMinVolumeFraction;
  arrayView2d< real64 const > m_phaseRelPermExponent;
  arrayView2d< real64 const > m_phaseRelPermMaxValue;
  arrayView1d< real64 const > m_volFracScale;
};

class BrooksCoreyRelativePermeability : public RelativePermeabilityBase
{
public:

  BrooksCoreyRelativePermeability( string const & name, dataRepository::Group * const parent );

//START_SPHINX_INCLUDE_00
  static string catalogName() { return "BrooksCoreyRelativePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrooksCoreyRelativePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

//START_SPHINX_INCLUDE_01
  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr char const * phaseMinVolumeFractionString() { return "phaseMinVolumeFraction"; }
    static constexpr char const * phaseRelPermExponentString() { return "phaseRelPermExponent"; }
    static constexpr char const * phaseRelPermMaxValueString() { return "phaseRelPermMaxValue"; }
    static constexpr char const * volFracScaleString() { return "volFracScale"; }
  } vieKeysBrooksCoreyRelativePermeability;
//END_SPHINX_INCLUDE_01

  arrayView2d< real64 const > getPhaseMinVolumeFraction() const override { return m_phaseMinVolumeFraction; };

protected:

  virtual void resizeFields( localIndex const size,
                             localIndex const numPts ) override;

  virtual void postInputInitialization() override;

//START_SPHINX_INCLUDE_02
  array2d< real64 > m_phaseMinVolumeFraction;
  array2d< real64 > m_phaseRelPermExponent;
  array2d< real64 > m_phaseRelPermMaxValue;

  array1d< real64 > m_volFracScale;
//END_SPHINX_INCLUDE_02
};

GEOS_HOST_DEVICE
inline void
BrooksCoreyRelativePermeabilityUpdate::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
           arraySlice1d< real64, constitutive::relperm::USD_PHASE - 2 > const & phaseTrappedVolFrac,
           arraySlice2d< real64, constitutive::relperm::USD_RELPERM - 2 > const & phaseRelPerm,
           arraySlice3d< real64, constitutive::relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseRelPerm_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  integer const numDir = 3;

  for( int dir = 0; dir < numDir; ++dir )
  {
    real64 const satScaleInv = 1.0 / m_volFracScale[dir];

    for( localIndex ip = 0; ip < numPhases(); ++ip )
    {
      real64 const satScaled = (phaseVolFraction[ip] - m_phaseMinVolumeFraction[dir][ip]) * satScaleInv;
      real64 const exponent  = m_phaseRelPermExponent[dir][ip];
      real64 const scale     = m_phaseRelPermMaxValue[dir][ip];

      if( satScaled > 0.0 && satScaled < 1.0 )
      {
        // intermediate value
        real64 const v = scale * pow( satScaled, exponent - 1.0 );

        phaseRelPerm[ip][dir] = v * satScaled;
        dPhaseRelPerm_dPhaseVolFrac[ip][ip][dir] = v * exponent * satScaleInv;
      }
      else
      {
        phaseRelPerm[ip][dir] = (satScaled <= 0.0) ? 0.0 : scale;
      }

      phaseTrappedVolFrac[ip] = LvArray::math::min( phaseVolFraction[ip], m_phaseMinVolumeFraction[dir][ip] );

    }
  }
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
