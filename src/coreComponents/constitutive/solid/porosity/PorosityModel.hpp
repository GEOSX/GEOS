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
 * @file PorosityModel.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_POROSITYMODEL_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_POROSITYMODEL_HPP_

#include "PorosityBase.hpp"

namespace geosx
{
namespace constitutive
{

class PorosityModelUpdates : public PorosityBaseUpdates
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

  PorosityModelUpdates( arrayView2d< real64 > const & newPorosity,
                        arrayView2d< real64 > const & oldPorosity,
                        arrayView2d< real64 > const & dPorosity_dPressure,
                        arrayView1d< real64 > const & referencePorosity,
                        arrayView2d< real64 > const & biotCoefficient,
                        real64 const & grainBulkModulus ):
    PorosityBaseUpdates( newPorosity,
                         oldPorosity,
                         dPorosity_dPressure,
                         referencePorosity ),
    m_biotCoefficient( grainBulkModulus ),
    m_grainBulkModulus( grainBulkModulus )
  {}

  /// Default copy constructor
  PorosityModelUpdates( PorosityModelUpdates const & ) = default;

  /// Default move constructor
  PorosityModelUpdates( PorosityModelUpdates && ) = default;

  /// Deleted copy assignment operator
  PorosityModelUpdates & operator=( POROSITYMODELUpdates const & ) = delete;

  /// Deleted move assignment operator
  PorosityModelUpdates & operator=( PorosityModelUpdates && ) = delete;


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void updatePorosity( localIndex const k,
                       localIndex const q,
                       real64 const & pressure,
                       real64 const & deltaPressure,
                       real64 const ( &strainIncrement )[6] ) const
  {
    real64 const biotSkeletonModulusInverse = ( m_biotCoefficient[k][q] - m_referencePorosity[k] ) / m_grainBulkModulus;

    porosity = m_oldPorosity[k][q] +
        + m_biotCoefficient[k][q] * LvArray::tensorOps::symTrace< 3 >( strainIncrement )
        + biotSkeletonModulusInverse * deltaPressure;

    savePorosity( k, q, porosity, biotSkeletonModulusInverse );

    real64 const biotTimesPressure = -biotCoefficient * ( pressure + deltaPressure );

    LvArray::tensorOps::symAddIdentity< 3 >( stress, biotTimesPressure );
  }

protected:
  arrayView2d< real64 > m_biotCoefficient;

  real64 m_grainBulkModulus;
};


class PorosityModel : public ConstitutiveBase
{
public:
  PorosityModel( string const & name, Group * const parent );

  virtual ~PorosityModel() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PorosityModel"; }

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


  using KernelWrapper = PorosityModelUpdates;

   /**
    * @brief Create an update kernel wrapper.
    * @return the wrapper
    */
   KernelWrapper createKernelUpdates()
   {
     return KernelWrapper( m_newPorosity,
                           m_oldPorosity,
                           m_dPorosity_dPressure,
                           m_referencePorosity,
                           m_biotCoefficient,
                           m_grainBulkModulus );
   }


protected:
  virtual void postProcessInput() override;

  array2d< real64 > m_newPorosity;

  array2d< real64 > m_oldPorosity;

  array2d< real64 > m_dPorosity_dPressure;

  array1d< real64 > m_referencePorosity;

  real64 m_defaultReferencePorosity;

  array2d< real64 > m_biotCoefficient;

  real64 m_grainBulkModulus;
};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_POROSITYMODEL_HPP_
