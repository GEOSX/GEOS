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
 * @file CoupledSolid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_COUPLEDSOLIDBASE_HPP_
#define GEOS_CONSTITUTIVE_SOLID_COUPLEDSOLIDBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/solid/SolidInternalEnergy.hpp"

namespace geos
{
namespace constitutive
{


class CoupledSolidBase : public ConstitutiveBase
{
public:

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  CoupledSolidBase( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~CoupledSolidBase() override;

  struct viewKeyStruct
  {
    static constexpr char const * solidModelNameString() { return "solidModelName"; }
    static constexpr char const * porosityModelNameString() { return "porosityModelName"; }
    static constexpr char const * permeabilityModelNameString() { return "permeabilityModelName"; }
    static constexpr char const * solidInternalEnergyModelNameString() { return "solidInternalEnergyModelName"; }
  };

  virtual std::vector< string > getSubRelationNames() const override final
  {
    std::vector< string > subRelationNames = { m_solidModelName,
                                               m_porosityModelName,
                                               m_permeabilityModelName };

    if( !m_solidInternalEnergyModelName.empty() )
    {
      subRelationNames.push_back( m_solidInternalEnergyModelName );
    }

    return subRelationNames;
  }

  /**
   * @brief get the old porosity.
   * return a constant arrayView2d to the old porosity
   */
  arrayView2d< real64 const > const getPorosity_n() const
  {
    return getBasePorosityModel().getPorosity_n();
  }

  /*
   * @brief get the new porosity.
   * return a constant arrayView2d to the new porosity
   */
  arrayView2d< real64 const > const getPorosity() const
  {
    return getBasePorosityModel().getPorosity();
  }

  /**
   * @brief get the reference porosity.
   * return a constant arrayView1d to the reference porosity
   */
  arrayView1d< real64 const > const getReferencePorosity() const
  {
    return getBasePorosityModel().getReferencePorosity();
  }

  /**
   * @brief get the dPorosity_dPressure.
   * return a constant arrayView2d to dPorosity_dPressure
   */
  arrayView2d< real64 const > const  getDporosity_dPressure() const
  { return getBasePorosityModel().dPorosity_dPressure(); }

  /**
   * @brief get the dPorosity_dTemperature.
   * return a constant arrayView2d to dPorosity_dTemperature
   */
  arrayView2d< real64 const > const  getDporosity_dTemperature() const
  { return getBasePorosityModel().dPorosity_dTemperature(); }


  /**
   * @brief get the old internal energy.
   * return a constant arrayView2d to the old internal energy
   */
  arrayView2d< real64 const > const getInternalEnergy_n() const
  {
    return getSolidInternalEnergyModel().getInternalEnergy_n();
  }

  /*
   * @brief get the internal energy.
   * return a constant arrayView2d to the internal energy
   */
  arrayView2d< real64 const > const getInternalEnergy() const
  {
    return getSolidInternalEnergyModel().getInternalEnergy();
  }

  /**
   * @brief get the dInternalEnergy_dTemperature.
   * return a constant arrayView2d to dInternalEnergy_dTemperature
   */
  arrayView2d< real64 const > const  getDinternalEnergy_dTemperature() const
  { return getSolidInternalEnergyModel().getDinternalEnergy_dTemperature(); }

  /*
   * @brief Utility function to scale the reference porosity (for instance, by net-to-gross)
   * @param[in] scalingFactors the vector of scaling factors for the reference porosity
   */
  void scaleReferencePorosity( arrayView1d< real64 const > scalingFactors ) const
  { getBasePorosityModel().scaleReferencePorosity( scalingFactors ); }

  /*
   * @brief get the current bulk modulus
   * return a constant arrayView1d to bulk modulus
   */
  arrayView1d< real64 const > const getBulkModulus() const
  {
    return getBaseSolidModel().getBulkModulus();
  }

  /*
   * @brief get the current bulk modulus
   * return a constant arrayView1d to bulk modulus
   */
  arrayView1d< real64 const > const getShearModulus() const
  {
    return getBaseSolidModel().getShearModulus();
  }

  /*
   * @brief get the current solid density
   * return a constant arrayView2d to solid density
   */
  arrayView2d< real64 const > const getDensity() const
  {
    return getBaseSolidModel().getDensity();
  }

  /*
   * @brief get the current solid effective stress
   * return a constant arrayView3d to effective stress in Voigt form
   */
  arrayView3d< real64 const, solid::STRESS_USD > const getStress() const
  {
    return getBaseSolidModel().getStress();
  }

  /*
   * @brief get the current biot coefficient
   * return a constant arrayView1d to biotCoefficient
   */
  arrayView1d< real64 const > const getBiotCoefficient() const
  {
    return getBasePorosityModel().getBiotCoefficient();
  }

  /**
   * @brief Const/non-mutable accessor for the mean total stress increment
   * (with respect to the previous converged time level) at the previous sequential iteration
   * @return Accessor
   */
  arrayView2d< real64 const > const getMeanTotalStressIncrement_k() const
  {
    return getBasePorosityModel().getMeanTotalStressIncrement_k();
  }

  /**
   * @brief Non-const accessor for the mean total stress increment at the previous sequential iteration
   * @return Accessor
   */
  arrayView1d< real64 > const getAverageMeanTotalStressIncrement_k()
  {
    return getBasePorosityModel().getAverageMeanTotalStressIncrement_k();
  }


  /**
   * @brief initialize the constitutive models fields.
   */
  virtual void initializeState() const
  {
    getBasePorosityModel().initializeState();
    getBasePermModel().initializeState();
  }

  virtual void saveConvergedState() const override final
  {
    getBasePorosityModel().saveConvergedState();
    if( !m_solidInternalEnergyModelName.empty() )
    {
      /// If the name is provided it has to be saved as well.
      getSolidInternalEnergyModel().saveConvergedState();
    }
  }

  /**
   * @brief ignore the porosity update (after initialization step)
   */
  virtual void ignoreConvergedState() const
  {
    getBasePorosityModel().ignoreConvergedState();
  }

  /**
   * @brief get a constant reference to the solid internal energy model
   * return a constant SolidInternalEnergy reference to the solid internal energy model
   */
  SolidInternalEnergy const & getSolidInternalEnergyModel() const
  { return this->getParent().template getGroup< SolidInternalEnergy >( m_solidInternalEnergyModelName ); }

  /**
   * @brief get a PorosityBase constant reference to the porosity model
   * return a constant PorosityBase reference to the porosity model
   */
  PorosityBase const & getBasePorosityModel() const
  { return this->getParent().template getGroup< PorosityBase >( m_porosityModelName ); }

  /**
   * @brief get a PorosityBase reference to the porosity model
   * return a PorosityBase reference to the porosity model
   */
  PorosityBase & getBasePorosityModel()
  { return this->getParent().template getGroup< PorosityBase >( m_porosityModelName ); }

protected:

  /// the name of the solid model
  string m_solidModelName;

  /// the name of the porosity model
  string m_porosityModelName;

  /// the name of the permeability model
  string m_permeabilityModelName;

  /// the name of the solid internal energy model
  string m_solidInternalEnergyModelName;

private:

  /**
   * @brief get a Permeability base constant reference to the permeability model
   * return a constant PermeabilityBase reference to the permeability model
   */
  PermeabilityBase const & getBasePermModel() const
  { return this->getParent().template getGroup< PermeabilityBase >( m_permeabilityModelName ); }

  /**
   *@brief get a SolidBase constant reference to the solid model
   * return a constant SolidBase reference to the solid model
   */
  SolidBase const & getBaseSolidModel() const
  { return this->getParent().template getGroup< SolidBase >( m_solidModelName ); }

};

}

}

#endif /* GEOS_CONSTITUTIVE_SOLID_COUPLEDSOLID_HPP_ */
