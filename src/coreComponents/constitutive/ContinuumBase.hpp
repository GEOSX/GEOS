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
 * @file ContinuumBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTINUUMBASE_HPP_
#define GEOS_CONSTITUTIVE_CONTINUUMBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Base class for all continuum constitutive kernel wrapper classes.
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
class ContinuumBaseUpdates
{
protected:
  /**
   * @brief constructor
   * @param[in] newStress The new stress data from the constitutive model class.
   * @param[in] oldStress The old stress data from the constitutive model class.
   * @param[in] wavespeed The wavespeed data from the constitutive model class.
   */
  ContinuumBaseUpdates( arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                        arrayView2d< real64 > const & density,
                        arrayView2d< real64 > const & wavespeed ):
    m_newStress( newStress ),
    m_oldStress( oldStress ),
    m_density( density ),
    m_wavespeed( wavespeed )
  {}

  /// Deleted default constructor
  ContinuumBaseUpdates() = delete;

  /**
   * @brief Copy Constructor
   * @param source Object to copy
   */
  ContinuumBaseUpdates( ContinuumBaseUpdates const & source ) = default;

  /**
   * @brief Move Constructor
   * @param source Object to move resources from
   */
  ContinuumBaseUpdates( ContinuumBaseUpdates && source ) = default;

  /// Deleted copy assignment operator
  ContinuumBaseUpdates & operator=( ContinuumBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  ContinuumBaseUpdates & operator=( ContinuumBaseUpdates && ) =  delete;

public:

  /// A reference the current material stress at quadrature points.
  arrayView3d< real64, solid::STRESS_USD > const m_newStress;

  /// A reference the previous material stress at quadrature points.
  arrayView3d< real64, solid::STRESS_USD > const m_oldStress;

  /// A reference the current material density at quadrature points.
  arrayView2d< real64 > m_density;

  /// A reference the current material wavespeed at quadrature points.
  arrayView2d< real64 > m_wavespeed;

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
};


/**
 * @class ContinuumBase
 * This class serves as the base class for solid constitutive models.
 */
class ContinuumBase : public constitutive::ConstitutiveBase
{
public:
  /**
   * @brief Constructor
   * @param name Name of the ContinuumBase object in the repository.
   * @param parent The parent group of the ContinuumBase object.
   */
  ContinuumBase( string const & name,
             Group * const parent );

  /**
   * Destructor
   */
  virtual ~ContinuumBase() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ContinuumBase";

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
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * stressString() { return "stress"; }                  ///< New stress key
    static constexpr char const * oldStressString() { return "oldStress"; }            ///< Old stress key
    static constexpr char const * densityString() { return "density"; }                ///< Density key
    static constexpr char const * defaultDensityString() { return "defaultDensity"; }  ///< Default density key
    static constexpr char const * wavespeedString() { return "wavespeed"; } 
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
   * @brief Number of elements storing solid data
   * @return Number of elements
   */
  localIndex numElem() const
  {
    return m_oldStress.size( 0 );
  }

  /**
   * @brief Number of quadrature points storing solid data
   * @return Number of quadrature points
   */
  localIndex numQuad() const
  {
    return m_oldStress.size( 1 );
  }

  /**
   * @name Accessors
   */
  ///@{

  /**
   * @brief Non-const/mutable accessor for stress
   * @return Accessor
   */
  arrayView3d< real64, solid::STRESS_USD > const getStress()
  {
    return m_newStress;
  }

  /**
   * @brief Const/non-mutable accessor for stress
   * @return Accessor
   */
  arrayView3d< real64 const, solid::STRESS_USD > const getStress() const
  {
    return m_newStress;
  }

  /**
   * @brief Non-const/mutable accessor for old stress
   * @return Accessor
   */
  arrayView3d< real64, solid::STRESS_USD > const getOldStress()
  {
    return m_oldStress;
  }

  /**
   * @brief Const/non-mutable accessor for old stress
   * @return Accessor
   */
  arrayView3d< real64 const, solid::STRESS_USD > const getOldStress() const
  {
    return m_oldStress;
  }

  /**
   * @brief Non-const/Mutable accessor for density.
   * @return Accessor
   */
  arrayView2d< real64 > const getDensity()
  {
    return m_density;
  }

  /**
   * @brief Const/non-mutable accessor for density
   * @return Accessor
   */
  arrayView2d< real64 const > const getDensity() const
  {
    return m_density;
  }

  /**
   * @brief Non-const/Mutable accessor for wavespeed.
   * @return Accessor
   */
  arrayView2d< real64 > const getWavespeed()
  {
    return m_wavespeed;
  }

  /**
   * @brief Const/non-mutable accessor for wavespeed.
   * @return Accessor
   */
  arrayView2d< real64 const > const getWavespeed() const
  {
    return m_wavespeed;
  }

  ///@}
  //

protected:

  /// Post-process XML input
  virtual void postInputInitialization() override;

  /// The current stress at a quadrature point (i.e. at timestep n, global newton iteration k)
  array3d< real64, solid::STRESS_PERMUTATION > m_newStress;

  /// The previous stress at a quadrature point (i.e. at timestep (n-1))
  array3d< real64, solid::STRESS_PERMUTATION > m_oldStress;

  /// The material density at a quadrature point.
  array2d< real64 > m_density;

  /// The default density for new allocations.
  real64 m_defaultDensity;

  /// The wavespeed at a quadrature point
  array2d< real64 > m_wavespeed;

  /// band-aid fix...going to have to remove this after we clean up
  /// initialization for constitutive models.
  bool m_postProcessed = false;
};

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_CONTINUUMBASE_HPP_ */
