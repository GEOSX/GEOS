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
 * @file PorousSolid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_
#define GEOS_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam SOLID_TYPE type of the porosity model
 */
template< typename SOLID_TYPE >
class PorousSolidUpdates : public CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >
{
public:

  using DiscretizationOps = typename SOLID_TYPE::KernelWrapper::DiscretizationOps;

  /**
   * @brief Constructor
   */
  PorousSolidUpdates( SOLID_TYPE const & solidModel,
                      BiotPorosity const & porosityModel,
                      ConstantPermeability const & permModel ):
    CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >( solidModel, porosityModel, permModel )
  {}

  GEOS_HOST_DEVICE
  virtual void updateStateFromPressureAndTemperature( localIndex const k,
                                                      localIndex const q,
                                                      real64 const & pressure,
                                                      real64 const & pressure_k,
                                                      real64 const & pressure_n,
                                                      real64 const & temperature,
                                                      real64 const & temperature_k,
                                                      real64 const & temperature_n ) const override final
  {
    updateBiotCoefficientAndAssignModuli( k );

    m_porosityUpdate.updateFixedStress( k, q,
                                        pressure, pressure_k, pressure_n,
                                        temperature, temperature_k, temperature_n );
  }

  GEOS_HOST_DEVICE
  void smallStrainUpdatePoromechanics( localIndex const k,
                                       localIndex const q,
                                       real64 const & timeIncrement,
                                       real64 const & pressure,
                                       real64 const & pressure_n,
                                       real64 const & temperature,
                                       real64 const & deltaTemperatureFromLastStep,
                                       real64 const ( &strainIncrement )[6],
                                       real64 ( & totalStress )[6],
                                       real64 ( & dTotalStress_dPressure )[6],
                                       real64 ( & dTotalStress_dTemperature )[6],
                                       DiscretizationOps & stiffness,
                                       integer const performStressInitialization,
                                       real64 & porosity,
                                       real64 & porosity_n,
                                       real64 & dPorosity_dVolStrain,
                                       real64 & dPorosity_dPressure,
                                       real64 & dPorosity_dTemperature,
                                       real64 & dSolidDensity_dPressure ) const
  {
    // Compute total stress increment and its derivative
    computeTotalStress( k,
                        q,
                        timeIncrement,
                        pressure,
                        temperature,
                        strainIncrement,
                        totalStress,
                        dTotalStress_dPressure,
                        dTotalStress_dTemperature,
                        stiffness );

    // Compute porosity and its derivatives
    real64 const deltaPressure = pressure - pressure_n;
    real64 porosityInit;
    computePorosity( k,
                     q,
                     deltaPressure,
                     deltaTemperatureFromLastStep,
                     strainIncrement,
                     porosity,
                     porosity_n,
                     porosityInit,
                     dPorosity_dVolStrain,
                     dPorosity_dPressure,
                     dPorosity_dTemperature );

    // skip porosity update when doing poromechanics initialization
    if( performStressInitialization )
    {
      porosity = porosityInit;
      dPorosity_dVolStrain = 0.0;
      dPorosity_dPressure = 0.0;
      dPorosity_dTemperature = 0.0;
    }

    // Save the derivative of solid density wrt pressure for the computation of the body force
    dSolidDensity_dPressure = m_porosityUpdate.dGrainDensity_dPressure( k );
  }

  GEOS_HOST_DEVICE
  void smallStrainUpdatePoromechanicsFixedStress( localIndex const k,
                                                  localIndex const q,
                                                  real64 const & timeIncrement,
                                                  real64 const & pressure,
                                                  real64 const & pressure_n,
                                                  real64 const & temperature,
                                                  real64 const & temperature_n,
                                                  real64 const ( &strainIncrement )[6],
                                                  real64 ( & totalStress )[6],
                                                  DiscretizationOps & stiffness ) const
  {
    real64 dTotalStress_dPressure[6]{};
    real64 dTotalStress_dTemperature[6]{};

    // Compute total stress increment and its derivative
    computeTotalStress( k,
                        q,
                        timeIncrement,
                        pressure,
                        temperature,
                        strainIncrement,
                        totalStress,
                        dTotalStress_dPressure, // To pass something here
                        dTotalStress_dTemperature, // To pass something here
                        stiffness );

    // Compute total stress increment for the porosity update
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );
    real64 const meanEffectiveStressIncrement = bulkModulus * ( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );
    real64 const biotCoefficient = m_porosityUpdate.getBiotCoefficient( k );
    real64 const thermalExpansionCoefficient = m_solidUpdate.getThermalExpansionCoefficient( k );
    real64 const meanTotalStressIncrement = meanEffectiveStressIncrement - biotCoefficient * ( pressure - pressure_n )
                                            - 3 * thermalExpansionCoefficient * bulkModulus * ( temperature - temperature_n );
    m_porosityUpdate.updateMeanTotalStressIncrement( k, q, meanTotalStressIncrement );
  }

  /**
   * @brief Return the stiffness at a given element (small-strain interface)
   *
   * @note If the material model has a strain-dependent material stiffness (e.g.
   * any plasticity, damage, or nonlinear elastic model) then this interface will
   * not work.  Users should instead use one of the interfaces where a strain
   * tensor is provided as input.
   *
   * @param k the element number
   * @param stiffness the stiffness array
   */
  GEOS_HOST_DEVICE
  inline
  void getElasticStiffness( localIndex const k, localIndex const q, real64 ( & stiffness )[6][6] ) const
  {
    m_solidUpdate.getElasticStiffness( k, q, stiffness );
  }

  /**
   * @brief Return the stiffness at a given element (small-strain interface)
   *
   * @param [in] k the element number
   * @param [out] biotCefficient the biot-coefficient
   */
  GEOS_HOST_DEVICE
  inline
  void getBiotCoefficient( localIndex const k, real64 & biotCoefficient ) const
  {
    biotCoefficient = m_porosityUpdate.getBiotCoefficient( k );
  }

  /**
   * @brief Return the stiffness at a given element (small-strain interface)
   *
   * @param [in] k the element number
   * @param [out] thermalExpansionCoefficient the thermal expansion coefficient
   */
  GEOS_HOST_DEVICE
  inline
  void getThermalExpansionCoefficient( localIndex const k, real64 & thermalExpansionCoefficient ) const
  {
    thermalExpansionCoefficient = m_solidUpdate.getThermalExpansionCoefficient( k );
  }

private:

  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_solidUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_porosityUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_permUpdate;


  GEOS_HOST_DEVICE
  inline
  void updateBiotCoefficientAndAssignModuli( localIndex const k ) const
  {
    // This call is not general like this.
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );
    real64 const shearModulus = m_solidUpdate.getShearModulus( k );

    m_porosityUpdate.updateBiotCoefficientAndAssignModuli( k, bulkModulus, shearModulus );
  }

  GEOS_HOST_DEVICE
  void updateThermalExpansionCoefficient( localIndex const k ) const
  {
    real64 const thermalExpansionCoefficient = m_solidUpdate.getThermalExpansionCoefficient( k );

    m_porosityUpdate.updateThermalExpansionCoefficient( k, thermalExpansionCoefficient );
  }

  GEOS_HOST_DEVICE
  void computePorosity( localIndex const k,
                        localIndex const q,
                        real64 const & deltaPressure,
                        real64 const & deltaTemperature,
                        real64 const ( &strainIncrement )[6],
                        real64 & porosity,
                        real64 & porosity_n,
                        real64 & porosityInit,
                        real64 & dPorosity_dVolStrain,
                        real64 & dPorosity_dPressure,
                        real64 & dPorosity_dTemperature ) const
  {
    m_porosityUpdate.updateFromPressureTemperatureAndStrain( k,
                                                             q,
                                                             deltaPressure,
                                                             deltaTemperature,
                                                             strainIncrement,
                                                             dPorosity_dVolStrain,
                                                             dPorosity_dPressure,
                                                             dPorosity_dTemperature );

    porosity = m_porosityUpdate.getPorosity( k, q );
    porosity_n = m_porosityUpdate.getPorosity_n( k, q );
    porosityInit = m_porosityUpdate.getInitialPorosity( k, q );
  }

  GEOS_HOST_DEVICE
  inline
  void computeTotalStress( localIndex const k,
                           localIndex const q,
                           real64 const & timeIncrement,
                           real64 const & pressure,
                           real64 const & temperature,
                           real64 const ( &strainIncrement )[6],
                           real64 ( & totalStress )[6],
                           real64 ( & dTotalStress_dPressure )[6],
                           real64 ( & dTotalStress_dTemperature )[6],
                           DiscretizationOps & stiffness ) const
  {
    updateBiotCoefficientAndAssignModuli( k );

    // Compute total stress increment and its derivative w.r.t. pressure
    m_solidUpdate.smallStrainUpdate( k,
                                     q,
                                     timeIncrement,
                                     strainIncrement,
                                     totalStress, // first effective stress increment accumulated
                                     stiffness );

    // Add the contributions of pressure and temperature to the total stress
    real64 const biotCoefficient = m_porosityUpdate.getBiotCoefficient( k );
    real64 const thermalExpansionCoefficient = m_solidUpdate.getThermalExpansionCoefficient( k );
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );
    real64 const thermalExpansionCoefficientTimesBulkModulus = thermalExpansionCoefficient * bulkModulus;

    LvArray::tensorOps::symAddIdentity< 3 >( totalStress, -biotCoefficient * pressure - 3 * thermalExpansionCoefficientTimesBulkModulus * temperature );

    // Compute derivatives of total stress
    dTotalStress_dPressure[0] = -biotCoefficient;
    dTotalStress_dPressure[1] = -biotCoefficient;
    dTotalStress_dPressure[2] = -biotCoefficient;
    dTotalStress_dPressure[3] = 0;
    dTotalStress_dPressure[4] = 0;
    dTotalStress_dPressure[5] = 0;

    dTotalStress_dTemperature[0] = -3 * thermalExpansionCoefficientTimesBulkModulus;
    dTotalStress_dTemperature[1] = -3 * thermalExpansionCoefficientTimesBulkModulus;
    dTotalStress_dTemperature[2] = -3 * thermalExpansionCoefficientTimesBulkModulus;
    dTotalStress_dTemperature[3] = 0;
    dTotalStress_dTemperature[4] = 0;
    dTotalStress_dTemperature[5] = 0;

  }
};

/**
 * @brief PorousSolidBase class used for dispatch of all Porous solids.
 */
class PorousSolidBase
{};

/**
 * @brief Class to represent a porous material for poromechanics simulations.
 * It is used as an interface to access all constitutive models relative to the properties of a porous material.
 *
 * @tparam SOLID_TYPE type of solid model
 */
template< typename SOLID_TYPE >
class PorousSolid : public CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >
{
public:

  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = PorousSolidUpdates< SOLID_TYPE >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  PorousSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~PorousSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Porous" ) + SOLID_TYPE::catalogName(); }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Create a instantiation of the PorousSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousSolidUpdates.
   */
  KernelWrapper createKernelUpdates() const
  {
    return KernelWrapper( getSolidModel(),
                          getPorosityModel(),
                          getPermModel() );
  }

  /**
   * @brief initialize the constitutive models fields.
   */
  virtual void initializeState() const override final;

  /**
   * @brief Const/non-mutable accessor for density
   * @return Accessor
   */
  arrayView2d< real64 const > const getDensity() const
  {
    return getSolidModel().getDensity();
  }

  /**
   * @brief Const/non-mutable accessor for biot coefficient
   * @return Accessor
   */
  arrayView1d< real64 const > const getBiotCoefficient() const
  {
    return getPorosityModel().getBiotCoefficient();
  }

  /**
   * @brief Const/non-mutable accessor for the mean stress increment at the previous sequential iteration
   * @return Accessor
   */
  arrayView2d< real64 const > const getMeanStressIncrement_k() const
  {
    return getPorosityModel().getMeanStressIncrement_k();
  }

  /**
   * @brief Non-const accessor for the mean stress increment at the previous sequential iteration
   * @return Accessor
   */
  arrayView1d< real64 > const getAverageMeanStressIncrement_k()
  {
    return getPorosityModel().getAverageMeanStressIncrement_k();
  }


private:
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getSolidModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPorosityModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPermModel;
};



}
} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_ */
