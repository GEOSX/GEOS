/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PorosityBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
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
  GEOSX_HOST_DEVICE
  localIndex numElems() const { return m_newPorosity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_newPorosity.size( 1 ); }

  PorosityBaseUpdates( arrayView2d< real64 > const & newPorosity,
                        arrayView2d< real64 > const & oldPorosity,
                        arrayView2d< real64 > const & dPorosity_dPressure,
                        arrayView1d< real64 > const & referencePorosity ):
    m_newPorosity( newPorosity ),
    m_oldPorosity( oldPorosity ),
    m_dPorosity_dPressure( dPorosity_dPressure ),
    m_referencePorosity ( referencePorosity )
  {}

  /// Default copy constructor
  PorosityBaseUpdates( PorosityBaseUpdates const & ) = default;

  /// Default move constructor
  PorosityBaseUpdates( PorosityBaseUpdates && ) = default;

  /// Deleted copy assignment operator
  PorosityBaseUpdates & operator=( POROSITYBASEUpdates const & ) = delete;

  /// Deleted move assignment operator
  PorosityBaseUpdates & operator=( PorosityBaseUpdates && ) = delete;


  /**
   * @brief Helper to save point stress back to m_newPorosity array
   *
   * This is mostly defined for improving code readability.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] porosity porosity to be saved to m_newPorosity[k][q]
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void savePorosity( localIndex const k,
                     localIndex const q,
                     real64 const & porosity,
                     real64 const & dPorosity_dPressure ) const
  {
    m_newPorosity[k][q] = porosity;
    m_dPorosity_dPressure[k][q] = dPorosity_dPressure;
  }

protected:
  arrayView2d< real64 > m_newPorosity;

  arrayView2d< real64 > m_oldPorosity;

  arrayView2d< real64 > m_dPorosity_dPressure;

  arrayView1d< real64 > m_referencePorosity;
};


class PorosityBase : public ConstitutiveBase
{
public:
  PorosityBase( string const & name, Group * const parent );

  virtual ~PorosityBase() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PorosityBase"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * newPorosityString() { return "porosity"; }
    static constexpr char const * oldPorosityString() { return "oldPorosity"; }
    static constexpr char const * dPorosity_dPressureString() { return "dPorosity_dPressure"; }
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }
    static constexpr char const * defaultRefererencePorosityString() { return "defaultReferencePorosity"; }
  } viewKeys;

  /**
   * @brief Const accessor for newPorosity.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getPorosity() const { return m_newPorosity; }

  /**
   * @brief Const/non-mutable accessor for oldPorosity.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getOldPorosity() const { return m_oldPorosity; }


  /**
   * @brief Non-Const/mutable accessor for oldPorosity
   * @return Accessor
   */
  arrayView2d< real64 > const getOldPorosity() { return m_oldPorosity; }


  /**
   * @brief Const/non-mutable accessor for dPorosity_dPressure
   * @return Accessor
   */
  arrayView2d< real64 const > const  dPorosity_dPressure() const { return m_dPorosity_dPressure; }


  using KernelWrapper = PorosityBaseUpdates;

   /**
    * @brief Create an update kernel wrapper.
    * @return the wrapper
    */
   KernelWrapper createKernelUpdates()
   {
     return KernelWrapper( m_newPorosity,
                           m_oldPorosity,
                           m_dPorosity_dPressure,
                           m_referencePorosity );
   }


protected:
  virtual void postProcessInput() override;

  array2d< real64 > m_newPorosity;

  array2d< real64 > m_oldPorosity;

  array2d< real64 > m_dPorosity_dPressure;

  array1d< real64 > m_referencePorosity;

  real64 m_defaultReferencePorosity;

};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_
