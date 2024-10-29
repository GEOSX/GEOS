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
 * @file SolidBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_

#include "constitutive/ContinuumBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Base class for all solid constitutive kernel wrapper classes.
 *
 * The responsibility of this base is to:
 *
 * 1) Contain views to state and parameter data for solid models.
 * 2) Specify an interface for state update functions.
 *
 * In general, the ArrayView data in the wrapper is specified to be of type
 * "arrayView<T> const" or "arrayView<T const> const". The "const-ness"
 * of the data indicates whether it is a parameter" or a state variable,
 * with the parameters being "T const" and state variables being "T".
 *
 * @note
 * If an allocation occurs on the underlying Array after a KernelWrapper is created,
 * then the ArrayView members of that KernelWrapper are silently invalid.
 */
class SolidBaseUpdates : public ContinuumBaseUpdates
{
protected:
  /**
   * @brief constructor
   * @param[in] newStress The new stress data from the constitutive model class.
   * @param[in] oldStress The old stress data from the constitutive model class.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] disableInelasticity Flag to disable inelastic response
   */
  SolidBaseUpdates( arrayView3d< real64, solid::STRESS_USD > const & newStress,
                    arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                    arrayView2d< real64 > const & density,
                    arrayView2d< real64 > const & wavespeed,
                    arrayView1d< real64 const > const & thermalExpansionCoefficient,
                    const bool & disableInelasticity ):
    ContinuumBaseUpdates( newStress,
                          oldStress,
                          density,
                          wavespeed ),
    m_thermalExpansionCoefficient( thermalExpansionCoefficient ),
    m_disableInelasticity( disableInelasticity )
  {}

  /// Deleted default constructor
  SolidBaseUpdates() = delete;

  /**
   * @brief Copy Constructor
   * @param source Object to copy
   */
  SolidBaseUpdates( SolidBaseUpdates const & source ) = default;

  /**
   * @brief Move Constructor
   * @param source Object to move resources from
   */
  SolidBaseUpdates( SolidBaseUpdates && source ) = default;

  /// Deleted copy assignment operator
  SolidBaseUpdates & operator=( SolidBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  SolidBaseUpdates & operator=( SolidBaseUpdates && ) =  delete;

public:

  /// A reference to the ArrayView holding the thermal expansion coefficient for each element.
  arrayView1d< real64 const > const m_thermalExpansionCoefficient;

  /// Flag to disable inelasticity
  const bool m_disableInelasticity;

  /**
   * @brief Get bulkModulus
   * @param[in] k Element index.
   * @return the bulkModulus of element k
   */
  GEOS_HOST_DEVICE
  virtual real64 getBulkModulus( localIndex const k ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_ERROR( "getBulkModulus() not implemented for this model" );

    return 0;
  }

  /**
   * @brief Get thermalExpansionCoefficient
   * @param[in] k Element index.
   * @return the thermalExpansionCoefficient of element k
   */
  GEOS_HOST_DEVICE
  real64 getThermalExpansionCoefficient( localIndex const k ) const
  {
    return m_thermalExpansionCoefficient[k];
  }

  /**
   * @brief Get shear modulus
   * @param[in] k Element index.
   * @return the shear modulus of element k
   */
  GEOS_HOST_DEVICE
  virtual real64 getShearModulus( localIndex const k ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_ERROR( "getShearModulus() not implemented for this model" );
    return 0;
  }

  /**
   * @name Update Interfaces: Stress and Stiffness
   *
   * We define a variety of interfaces for constitutive models using different
   * strain theories.  Derived classes only need to implement the subset of interfaces
   * most relevant to them.
   *
   * This group of interfaces returns stress and stiffness simultaneously, and
   * are most useful for implicit finite element formulations.
   */
  ///@{

  /**
   * @brief Small strain update.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] timeIncrement time increment for rate-dependent models.
   * @param[in] strainIncrement Strain increment in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New tangent stiffness value
   */
  GEOS_HOST_DEVICE
  /**
   * this function is not virtual to avoid a compilation warning with nvcc.
   */
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( & stress )[6],
                          real64 ( & stiffness )[6][6] ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( timeIncrement );
    GEOS_UNUSED_VAR( strainIncrement );
    GEOS_UNUSED_VAR( stress );
    GEOS_UNUSED_VAR( stiffness );
    GEOS_ERROR( "smallStrainUpdate() not implemented for this model" );
  }

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_ElasticOnly( localIndex const k,
                                              localIndex const q,
                                              real64 const & timeIncrement,
                                              real64 const ( &strainIncrement )[6],
                                              real64 ( & stress )[6],
                                              real64 ( & stiffness )[6][6] ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( timeIncrement );
    GEOS_UNUSED_VAR( strainIncrement );
    GEOS_UNUSED_VAR( stress );
    GEOS_UNUSED_VAR( stiffness );
    GEOS_ERROR( "smallStrainUpdate_ElasticOnly() not implemented for this model, or the model is already elastic." );
  }

  /**
   * @brief Small strain, stateless update.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] totalStrain Total strain in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   * @param[out] stiffness New tangent stiffness value
   */
  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( & stress )[6],
                                         real64 ( & stiffness )[6][6] ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( totalStrain );
    GEOS_UNUSED_VAR( stress );
    GEOS_UNUSED_VAR( stiffness );
    GEOS_ERROR( "smallStrainNoStateUpdate() not implemented for this model" );
  }

  ///@}
  /**
   * @name Update Interfaces: Stress-Only
   *
   * We define a variety of interfaces for constitutive models using different
   * strain theories.  Derived classes only need to implement the subset of interfaces
   * most relevant to them.
   *
   * This group of interfaces returns stress only, with no stiffness, and
   * are most useful for explicit finite element formulations.
   *
   * @note
   * The base class versions implement a naive call to the versions
   * of the updates that returns both stress and stiffness, and then discards
   * the stiffness.  Derived classes can implement an optimized version to
   * avoid extranenous work, but we delegate this detail to them.
   */
  ///@{

  /**
   * @brief Small strain update, returning only stress.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] timeIncrement time increment for rate-dependent models.
   * @param[in] strainIncrement Strain increment in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   */
  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( & strainIncrement )[6],
                                             real64 ( & stress )[6] ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( timeIncrement );
    GEOS_UNUSED_VAR( strainIncrement );
    GEOS_UNUSED_VAR( stress );
    GEOS_ERROR( "smallStrainUpdate_StressOnly() not implemented for this model" );
  }


  /**
   * @brief Small strain update overload with rotations, returning only stress for material point method solver.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] timeIncrement time increment for rate-dependent models.
   * @param[in] beginningRotation rotation matrix at beginning of time step (used to unrotate state variables in certain models)
   * @param[in] endRotation rotation matrix at end of time step
   * @param[in] strainIncrement Strain increment in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   */
  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( & beginningRotation )[3][3],
                                             real64 const ( & endRotation )[3][3],
                                             real64 const ( & strainIncrement )[6],
                                             real64 ( & stress )[6] ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( timeIncrement );
    GEOS_UNUSED_VAR( beginningRotation );
    GEOS_UNUSED_VAR( endRotation );
    GEOS_UNUSED_VAR( strainIncrement );
    GEOS_UNUSED_VAR( stress );
    GEOS_ERROR( "smallStrainUpdate_StressOnly() not implemented for this model" );
  }


  /**
   * @brief Small strain, stateless update, returning only stress.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] totalStrain total strain in Voight notation (linearized strain)
   * @param[out] stress New stress value (Cauchy stress)
   */
  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( & stress )[6] ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( totalStrain );
    GEOS_UNUSED_VAR( stress );
    GEOS_ERROR( "smallStrainNoStateUpdate_StressOnly() not implemented for this model" );
  }

  /**
   * @brief Helper to save point stress back to m_newStress array
   *
   * This is mostly defined for improving code readability.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] stress Stress to be save to m_newStress[k][q]
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void saveStress( localIndex const k,
                   localIndex const q,
                   real64 const ( &stress )[6] ) const
  {
    LvArray::tensorOps::copy< 6 >( m_newStress[k][q], stress );
  }

  /**
   * @brief Save converged state data at index (k,q)
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   */
  GEOS_HOST_DEVICE
  inline
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const
  {
    LvArray::tensorOps::copy< 6 >( m_oldStress[k][q], m_newStress[k][q] );
  }

  /**
   * @brief Return the current elastic strain at a given material point (small-strain interface)
   *
   * @param k the element inex
   * @param q the quadrature index
   * @param elasticStrain Current elastic strain
   */
  GEOS_HOST_DEVICE
  virtual void getElasticStrain( localIndex const k,
                                 localIndex const q,
                                 real64 ( & elasticStrain )[6] ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( elasticStrain );
    GEOS_ERROR( "getElasticStrain() not implemented for this model" );
  }

  /**
   * @brief Perform a viscous (rate-dependent) state update
   *
   * @param beta time-dependent parameter
   */
  GEOS_HOST_DEVICE
  virtual void viscousStateUpdate( localIndex const k,
                                   localIndex const q,
                                   real64 beta ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( beta );
    GEOS_ERROR( "viscousStateUpdate() not implemented for this model" );
  }

  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6] ) const 
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( FminusI );
    GEOS_UNUSED_VAR( stress );
    GEOS_ERROR( "hyperUpdate() not implemented for this model" );
  }

  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6],
                            real64 ( & stiffness )[6][6] ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( FminusI );
    GEOS_UNUSED_VAR( stress );
    GEOS_UNUSED_VAR( stiffness );
    GEOS_ERROR( "hyperUpdate() not implemented for this model" );
  }                     

  /**
   * @brief Return the strain energy density at a given material point
   *
   * @param k the element inex
   * @param q the quadrature index
   * @return Strain energy density
   */
  GEOS_HOST_DEVICE
  virtual real64 getStrainEnergyDensity( localIndex const k,
                                         localIndex const q ) const
  {
    auto const & stress = m_newStress[k][q];

    real64 strain[6]{};
    getElasticStrain( k, q, strain );

    real64 energy = 0;

    for( localIndex i=0; i<6; ++i )
    {
      energy += stress[i]*strain[i];  // contraction sigma:epsilon
    }
    energy *= 0.5;

    GEOS_ASSERT_MSG( energy >= 0.0, "negative strain energy density detected" );

    return energy;
  }

  /**
   * @brief Return the stiffness at a given element (small-strain interface)
   *
   * @note If the material model has a strain-dependent material stiffness (e.g.
   * any plasticity, damage, or nonlinear elastic model) then this interface will
   * not work.  Users should instead use one of the interfaces where a strain
   * tensor is provided as input.
   *
   * @note Given the limitations above, this function may be removed from the
   * public interface in the future.  Direct use in physics
   * solvers is discouraged.
   *
   * @param k the element number
   * @param stiffness the stiffness array
   */
  GEOS_HOST_DEVICE
  virtual void getElasticStiffness( localIndex const k, localIndex const q, real64 ( & stiffness )[6][6] ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( stiffness );
    GEOS_ERROR( "getElasticStiffness() not implemented for this model" );
  }

  /**
   * @brief Perform a finite-difference stiffness computation
   *
   * This method uses stress evaluations and finite differencing to
   * approximate the 6x6 stiffness matrix.
   *
   * @note This method only works for models providing the smallStrainUpdate
   * method returning a 6x6 stiffness, as it will primarily be used to check
   * the hand coded tangent against a finite difference reference.
   * A similar method would need to be implemented to check compressed stiffness,
   * stress-only, or finite-strain interfaces.
   *
   * @param k the element number
   * @param q the quadrature index
   * @param strainIncrement strain increment (on top of which a FD perturbation will be added)
   * @param stiffnessFD finite different stiffness approximation
   */
  GEOS_HOST_DEVICE
  void computeSmallStrainFiniteDifferenceStiffness( localIndex k,
                                                    localIndex q,
                                                    real64 const & timeIncrement,
                                                    real64 const ( &strainIncrement )[6],
                                                    real64 ( & stiffnessFD )[6][6] ) const
  {
    real64 stiffness[6][6]{};      // coded stiffness
    real64 stress[6]{};            // original stress
    real64 stressFD[6]{};          // perturbed stress
    real64 strainIncrementFD[6]{}; // perturbed strain
    real64 norm = 0;               // norm for scaling (note: method is fragile w.r.t. scaling)

    for( localIndex i=0; i<6; ++i )
    {
      strainIncrementFD[i] = strainIncrement[i];
      norm += fabs( strainIncrement[i] );
    }

    real64 eps = 1e-4*norm;     // finite difference perturbation

    smallStrainUpdate( k,
                       q, 
                       timeIncrement,
                       strainIncrement, 
                       stress, 
                       stiffness );

    for( localIndex i=0; i<6; ++i )
    {
      strainIncrementFD[i] += eps;

      if( i>0 )
      {
        strainIncrementFD[i-1] -= eps;
      }

      smallStrainUpdate( k, 
                         q, 
                         timeIncrement,
                         strainIncrementFD, 
                         stressFD, 
                         stiffnessFD );

      for( localIndex j=0; j<6; ++j )
      {
        stiffnessFD[j][i] = (stressFD[j]-stress[j])/eps;
      }
    }

    return;
  }

  /**
   * @brief Perform a finite-difference check of the stiffness computation
   *
   * This method uses several stress evaluations and finite differencing to
   * approximate the 6x6 stiffness matrix, and then computes an error between
   * the coded stiffness method and the finite difference version.
   *
   * @note This method only works for models providing the smallStrainUpdate
   * method returning a 6x6 stiffness.
   *
   * @param k the element number
   * @param q the quadrature index
   * @param strainIncrement strain increment (on top of which a FD perturbation will be added)
   */
  GEOS_HOST_DEVICE
  bool checkSmallStrainStiffness( localIndex k,
                                  localIndex q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  bool print = false ) const
  {
    real64 stiffness[6][6]{};     // coded stiffness
    real64 stiffnessFD[6][6]{};   // finite difference approximation
    real64 stress[6]{};           // original stress

    smallStrainUpdate( k, 
                       q, 
                       timeIncrement,
                       strainIncrement, 
                       stress, 
                       stiffness );
    computeSmallStrainFiniteDifferenceStiffness( k, 
                                                 q, 
                                                 timeIncrement, 
                                                 strainIncrement, 
                                                 stiffnessFD );

    // compute relative error between two versions

    real64 error = 0;
    real64 norm = 0;

    for( localIndex i=0; i<6; ++i )
    {
      for( localIndex j=0; j<6; ++j )
      {
        error += fabs( stiffnessFD[i][j]-stiffness[i][j] );
        norm += fabs( stiffnessFD[i][j] );
      }
    }
    error /= norm;

    // optional printing for debugging purposes

    if( print )
    {
      for( localIndex i=0; i<6; ++i )
      {
        for( localIndex j=0; j<6; ++j )
        {
          printf( "[%8.1e vs %8.1e] ", stiffnessFD[i][j], stiffness[i][j] );
        }
        printf( "\n" );
      }
    }

    return (error < 1e-3);
  }

};


/**
 * @class SolidBase
 * This class serves as the base class for solid constitutive models.
 */
class SolidBase : public constitutive::ContinuumBase
{
public:
  /**
   * @brief Constructor
   * @param name Name of the SolidBase object in the repository.
   * @param parent The parent group of the SolidBase object.
   */
  SolidBase( string const & name,
             Group * const parent );

  /**
   * Destructor
   */
  virtual ~SolidBase() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "SolidBase";

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// Keys for data in this class
  struct viewKeyStruct : public ContinuumBase::viewKeyStruct
  {
    static constexpr char const * thermalExpansionCoefficientString() { return "thermalExpansionCoefficient"; } // Thermal expansion
                                                                                                                // coefficient key
    static constexpr char const * defaultThermalExpansionCoefficientString() { return "defaultDrainedLinearTEC"; } // Default
                                                                                                                   // drained
                                                                                                                   // linear
                                                                                                                   // thermal
                                                                                                                   // expansion
                                                                                                                   // coefficient
                                                                                                                   // key
  };

  /**
   * @brief Allocate constitutive arrays
   * @param parent Object's parent group (element subregion)
   * @param numConstitutivePointsPerParentIndex Number of quadrature points per element
   */
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// Save state data in preparation for next timestep
  virtual void saveConvergedState() const override;

  /**
   * @brief Enable/disable inelasticity
   * @param flag Flag to disable (if true) or enable (if false) inelastic response
   */
  void disableInelasticity( bool const flag )
  {
    m_disableInelasticity = flag;
  }

  /**
   *@brief Get bulkModulus vector
   *@return the bulkModulus of all elements
   */
  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getBulkModulus( ) const
  {
    GEOS_ERROR( "getBulkModulus() not implemented for this model" );

    array1d< real64 > out;
    return out.toViewConst();
  }

  /**
   *@brief Get shearModulus vector
   *@return the shearModulus of all elements
   */
  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getShearModulus( ) const
  {
    GEOS_ERROR( "getShearModulus() not implemented for this model" );

    array1d< real64 > out;
    return out.toViewConst();
  }

protected:

  /// Post-process XML input
  virtual void postInputInitialization() override;

  /// The thermal expansion coefficient for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_thermalExpansionCoefficient;

  /// The default value of the thermal expansion coefficient for any new allocations.
  real64 m_defaultThermalExpansionCoefficient = 0;

  /// Flag to disable inelasticity (plasticity, damage, etc.)
  bool m_disableInelasticity = false;
};

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDBASE_HPP_ */
