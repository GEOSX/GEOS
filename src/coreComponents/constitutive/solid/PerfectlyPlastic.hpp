/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file PerfectlyPlastic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_PERFECTLYPLASTIC_HPP
#define GEOSX_CONSTITUTIVE_SOLID_PERFECTLYPLASTIC_HPP

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class PerfectlyPlasticUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class PerfectlyPlasticUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] yieldStress The ArrayView holding the yield stress for each element.
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  PerfectlyPlasticUpdates( arrayView1d< real64 const > const & yieldStress,
                           arrayView1d< real64 const > const & bulkModulus,
                           arrayView1d< real64 const > const & shearModulus,
                           arrayView1d< real64 const > const & thermalExpansionCoefficient,
                           arrayView3d< real64, solid::STRESS_USD > const & newStress,
                           arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                           bool const & disableInelasticity ):
    ElasticIsotropicUpdates( bulkModulus, shearModulus, thermalExpansionCoefficient, newStress, oldStress, disableInelasticity ),
    m_yieldStress( yieldStress )
  {}

  /// Default copy constructor
  PerfectlyPlasticUpdates( PerfectlyPlasticUpdates const & ) = default;

  /// Default move constructor
  PerfectlyPlasticUpdates( PerfectlyPlasticUpdates && ) = default;

  /// Deleted default constructor
  PerfectlyPlasticUpdates() = delete;

  /// Deleted copy assignment operator
  PerfectlyPlasticUpdates & operator=( PerfectlyPlasticUpdates const & ) = delete;

  /// Deleted move assignment operator
  PerfectlyPlasticUpdates & operator=( PerfectlyPlasticUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using ElasticIsotropicUpdates::smallStrainUpdate;

  GEOS_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( &stress )[6],
                          real64 ( &stiffness )[6][6] ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
  }

private:
  /// A reference to the ArrayView holding the yield stress for each element.
  arrayView1d< real64 const > const m_yieldStress;

};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void PerfectlyPlasticUpdates::smallStrainUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const & timeIncrement,
                                                 real64 const ( &strainIncrement )[6],
                                                 real64 ( & stress )[6],
                                                 real64 ( & stiffness )[6][6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );

  if( m_disableInelasticity )
  {
    return;
  }

  // decompose into mean (P) and von Mises (Q) stress invariants
  real64 trialP;
  real64 trialQ;
  real64 deviator[6];
  twoInvariant::stressDecomposition( stress,
                                     trialP,
                                     trialQ,
                                     deviator );

  // check yield function
  real64 yield = trialQ / m_yieldStress[k];
  if( yield < 1.0 ) // elasticity
  {
    return;
  }

  // else, plasticity
  // re-construct stress = P*eye + sqrt(2/3)*Q*nhat
  twoInvariant::stressRecomposition( trialP,
                                     m_yieldStress[k],
                                     deviator,
                                     stress );

  // TODO: construct consistent tangent stiffness (assuming additional yielding will occur?)
  // We don't modify the stiffness for now...

  // save new stress and return
  saveStress( k, q, stress );
  return;
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void PerfectlyPlasticUpdates::smallStrainUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const & timeIncrement,
                                                 real64 const ( &strainIncrement )[6],
                                                 real64 ( & stress )[6],
                                                 DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness.m_c );
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE // TODO: Is there a way to not have to re-write the constitutive model twice for regular and StressOnly?
void PerfectlyPlasticUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                            localIndex const q,
                                                            real64 const & timeIncrement,
                                                            real64 const ( &strainIncrement )[6],
                                                            real64 ( & stress )[6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate_StressOnly( k, q, timeIncrement, strainIncrement, stress );

  if( m_disableInelasticity )
  {
    return;
  }

  // decompose into mean (P) and von Mises (Q) stress invariants
  real64 trialP;
  real64 trialQ;
  real64 deviator[6];
  twoInvariant::stressDecomposition( stress,
                                     trialP,
                                     trialQ,
                                     deviator );

  // check yield function
  real64 yield = trialQ / m_yieldStress[k];
  if( yield < 1.0 ) // elasticity
  {
    return;
  }

  // else, plasticity
  // re-construct stress = P*eye + sqrt(2/3)*Q*nhat
  twoInvariant::stressRecomposition( trialP,
                                     m_yieldStress[k],
                                     deviator,
                                     stress );

  // save new stress and return
  saveStress( k, q, stress );
  return;
}



/**
 * @class PerfectlyPlastic
 *
 * Perfectly plastic material model.
 */
class PerfectlyPlastic : public ElasticIsotropic
{
public:

  /// @typedef Alias for PerfectlyPlasticUpdates
  using KernelWrapper = PerfectlyPlasticUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  PerfectlyPlastic( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~PerfectlyPlastic() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "PerfectlyPlastic";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    /// string/key for default yield stress
    static constexpr char const * defaultYieldStressString() { return "defaultYieldStress"; }

    /// string/key for elemental yield stress
    static constexpr char const * yieldStressString() { return "yieldStress"; }
  };

  /**
   * @brief Create a instantiation of the PerfectlyPlasticUpdate class that refers to the data in this.
   * @return An instantiation of PerfectlyPlasticUpdate.
   */
  PerfectlyPlasticUpdates createKernelUpdates() const
  {
    return PerfectlyPlasticUpdates( m_yieldStress,
                                    m_bulkModulus,
                                    m_shearModulus,
                                    m_thermalExpansionCoefficient,
                                    m_newStress,
                                    m_oldStress,
                                    m_disableInelasticity );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_yieldStress,
                          m_bulkModulus,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_disableInelasticity );
  }


protected:
  virtual void postInputInitialization() override;

  /// Material parameter: The default value of yield stress
  real64 m_defaultYieldStress;

  /// Material parameter: The yield stress for each element
  array1d< real64 > m_yieldStress;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_PERFECTLYPLASTIC_HPP_ */
