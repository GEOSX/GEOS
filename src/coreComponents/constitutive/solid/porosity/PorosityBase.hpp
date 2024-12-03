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
 * @file PorosityBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_
#define GEOS_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{
namespace constitutive
{

class PorosityBaseUpdates
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_newPorosity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_newPorosity.size( 1 ); }

  PorosityBaseUpdates( arrayView2d< real64 > const & newPorosity,
                       arrayView2d< real64 const > const & porosity_n,
                       arrayView2d< real64 > const & dPorosity_dPressure,
                       arrayView2d< real64 > const & dPorosity_dTemperature,
                       arrayView2d< real64 const > const & initialPorosity,
                       arrayView1d< real64 const > const & referencePorosity ):
    m_newPorosity( newPorosity ),
    m_porosity_n( porosity_n ),
    m_dPorosity_dPressure( dPorosity_dPressure ),
    m_dPorosity_dTemperature( dPorosity_dTemperature ),
    m_initialPorosity( initialPorosity ),
    m_referencePorosity ( referencePorosity )
  {}

  /**
   * @brief Helper to save porosity back to m_newPorosity array
   *
   * This is mostly defined for improving code readability.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] porosity porosity to be saved to m_newPorosity[k][q]
   * @param[in] dPorosity_dPressure porosity derivative w.r.t pressure to be saved to m_dPorosity_dPressure[k][q]
   * @param[in] dPorosity_dTemperature porosity derivative w.r.t temperature to be saved to m_dPorosity_dTemperature[k][q]
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void savePorosity( localIndex const k,
                     localIndex const q,
                     real64 const & porosity,
                     real64 const & dPorosity_dPressure,
                     real64 const & dPorosity_dTemperature ) const
  {
    m_newPorosity[k][q] = porosity;
    m_dPorosity_dPressure[k][q] = dPorosity_dPressure;
    m_dPorosity_dTemperature[k][q] = dPorosity_dTemperature;
  }

  GEOS_HOST_DEVICE
  inline
  real64 getPorosity( localIndex const k,
                      localIndex const q ) const
  {
    return m_newPorosity[k][q];
  }


  GEOS_HOST_DEVICE
  inline
  real64 getPorosity_n( localIndex const k,
                        localIndex const q ) const
  {
    return m_porosity_n[k][q];
  }

  GEOS_HOST_DEVICE
  inline
  real64 getInitialPorosity( localIndex const k,
                             localIndex const q ) const
  {
    return m_initialPorosity[k][q];
  }

protected:

  /// New value of porosity
  arrayView2d< real64 > const m_newPorosity;

  /// Value of porosity at the previous time step
  arrayView2d< real64 const > const m_porosity_n;

  /// Derivative of porosity wrt pressure
  arrayView2d< real64 > const m_dPorosity_dPressure;

  /// Derivative of porosity wrt temperature
  arrayView2d< real64 > const m_dPorosity_dTemperature;

  /// Initial porosity
  arrayView2d< real64 const > const m_initialPorosity;

  /// Reference porosity
  arrayView1d< real64 const > const m_referencePorosity;
};


class PorosityBase : public ConstitutiveBase
{
public:
  PorosityBase( string const & name, Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * defaultReferencePorosityString() { return "defaultReferencePorosity"; }
  } viewKeys;

  /**
   * @brief Number of elements storing solid data
   * @return Number of elements
   */
  localIndex numElem() const
  {
    return m_porosity_n.size( 0 );
  }

  /**
   * @brief Number of quadrature points storing solid data
   * @return Number of quadrature points
   */
  localIndex numQuad() const
  {
    return m_porosity_n.size( 1 );
  }


  /**
   * @brief Const accessor for newPorosity.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getPorosity() const { return m_newPorosity; }

  /**
   * @brief Const/non-mutable accessor for porosity_n.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getPorosity_n() const { return m_porosity_n; }


  /**
   * @brief Non-Const/mutable accessor for porosity_n
   * @return Accessor
   */
  arrayView2d< real64 > const getPorosity_n() { return m_porosity_n; }


  /**
   * @brief Const/non-mutable accessor for referencePorosity.
   * @return Accessor
   */
  arrayView1d< real64 const > const  getReferencePorosity() const { return m_referencePorosity; }


  /**
   * @brief Const/non-mutable accessor for dPorosity_dPressure
   * @return Accessor
   */
  arrayView2d< real64 const > const  dPorosity_dPressure() const { return m_dPorosity_dPressure; }

  /**
   * @brief Const/non-mutable accessor for dPorosity_dTemperature
   * @return Accessor
   */
  arrayView2d< real64 const > const  dPorosity_dTemperature() const { return m_dPorosity_dTemperature; }

  /**
   * @brief Utility function to scale the reference porosity (for instance, by net-to-gross)
   * @param[in] scalingFactors the vector of scaling factors (one value per cell) for the reference porosity
   */
  void scaleReferencePorosity( arrayView1d< real64 const > scalingFactors ) const;

  /// Save state data in preparation for next timestep
  virtual void saveConvergedState() const override;

  /// Ignore the porosity update and return to the state of the system
  /// This is useful after the initialization step
  virtual void ignoreConvergedState() const;

  /**
   * @brief Initialize newPorosity and porosity_n.
   */
  virtual void initializeState() const;

  virtual arrayView1d< real64 const > const getBiotCoefficient() const
  {
    GEOS_ERROR( "getBiotCoefficient() not implemented for this model" );

    array1d< real64 > out;
    return out.toViewConst();
  }

  /**
   * @brief Const/non-mutable accessor for the mean total stress increment at the previous sequential iteration
   * @return Accessor
   */
  virtual arrayView2d< real64 const > const getMeanTotalStressIncrement_k() const
  {
    GEOS_ERROR( "getMeanTotalStressIncrement_k() not implemented for this model" );

    array2d< real64 > out;
    return out.toViewConst();
  }

  /**
   * @brief Non-const accessor for the average mean total stress increment at the previous sequential iteration
   * @return Accessor
   */
  virtual arrayView1d< real64 > const getAverageMeanTotalStressIncrement_k()
  {
    GEOS_ERROR( "getAverageMeanTotalStressIncrement_k() not implemented for this model" );

    array1d< real64 > out;
    return out.toView();
  }

  using KernelWrapper = PorosityBaseUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates()
  {
    return KernelWrapper( m_newPorosity,
                          m_porosity_n,
                          m_dPorosity_dPressure,
                          m_dPorosity_dTemperature,
                          m_initialPorosity,
                          m_referencePorosity );
  }


protected:
  virtual void postInputInitialization() override;

  array2d< real64 > m_newPorosity;

  array2d< real64 > m_porosity_n;

  array2d< real64 > m_dPorosity_dPressure;

  array2d< real64 > m_dPorosity_dTemperature;

  array2d< real64 > m_initialPorosity;

  array1d< real64 > m_referencePorosity;

  real64 m_defaultReferencePorosity;

};

} /* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_
