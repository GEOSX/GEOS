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
 *  @file DelftEgg.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_DELFTEGG_HPP
#define GEOS_CONSTITUTIVE_SOLID_DELFTEGG_HPP

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
 * @class DelftEggUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class DelftEggUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] recompressionIndex          The ArrayView holding the recompression index data for each element.
   * @param[in] virginCompressionIndex      The ArrayView holding the virgin compression index data for each element.
   * @param[in] cslSlope                    The ArrayView holding the slope of the critical state line data for each element.
   * @param[in] shapeParameter              The ArrayView holding the shape parameter of yield surface data for each element.
   * @param[in] newPreConsolidationPressure The ArrayView holding the new preconsolidation pressure data for each qudrature point.
   * @param[in] oldPreConsolidationPressure The ArrayView holding old preconsolidation pressure data from the previous converged state for
   * each quadrature point.
   * @param[in] bulkModulus                 The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus                The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress                   The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress                   The ArrayView holding the old stress data from the previous converged state for each point
   * @param[in] disableInelasticity         Flag to disable plastic response
   */
  DelftEggUpdates( arrayView1d< real64 const > const & recompressionIndex,
                   arrayView1d< real64 const > const & virginCompressionIndex,
                   arrayView1d< real64 const > const & cslSlope,
                   arrayView1d< real64 const > const & shapeParameter,
                   arrayView2d< real64 > const & newPreConsolidationPressure,
                   arrayView2d< real64 > const & oldPreConsolidationPressure,
                   arrayView1d< real64 const > const & bulkModulus,
                   arrayView1d< real64 const > const & shearModulus,
                   arrayView1d< real64 const > const & thermalExpansionCoefficient,
                   arrayView3d< real64, solid::STRESS_USD > const & newStress,
                   arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                   const bool & disableInelasticity ):
    ElasticIsotropicUpdates( bulkModulus, shearModulus, thermalExpansionCoefficient, newStress, oldStress, disableInelasticity ),
    m_recompressionIndex( recompressionIndex ),
    m_virginCompressionIndex( virginCompressionIndex ),
    m_cslSlope( cslSlope ),
    m_shapeParameter( shapeParameter ),
    m_newPreConsolidationPressure( newPreConsolidationPressure ),
    m_oldPreConsolidationPressure( oldPreConsolidationPressure )
  {}

  /// Default copy constructor
  DelftEggUpdates( DelftEggUpdates const & ) = default;

  /// Default move constructor
  DelftEggUpdates( DelftEggUpdates && ) = default;

  /// Deleted default constructor
  DelftEggUpdates() = delete;

  /// Deleted copy assignment operator
  DelftEggUpdates & operator=( DelftEggUpdates const & ) = delete;

  /// Deleted move assignment operator
  DelftEggUpdates & operator=( DelftEggUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using ElasticIsotropicUpdates::smallStrainUpdate;

  GEOS_HOST_DEVICE
  void evaluateYield( real64 const p,
                      real64 const q,
                      real64 const pc,
                      real64 const M,
                      real64 const alpha,
                      real64 const Cc,
                      real64 const Cr,
                      real64 const bulkModulus,
                      real64 const mu,
                      real64 & f,
                      real64 & df_dp,
                      real64 & df_dq,
                      real64 & df_dpc,
                      real64 & df_dp_dve,
                      real64 & df_dq_dse ) const;



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
                                  DiscretizationOps & stiffness ) const final;

  GEOS_HOST_DEVICE
  inline
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
    m_oldPreConsolidationPressure[k][q] = m_newPreConsolidationPressure[k][q];
  }

private:

  /// A reference to the ArrayView holding the recompression index for each element.
  arrayView1d< real64 const > const m_recompressionIndex;

  /// A reference to the ArrayView holding the virgin compression index for each element.
  arrayView1d< real64 const > const m_virginCompressionIndex;

  /// A reference to the ArrayView holding the slope of the critical state line for each element.
  arrayView1d< real64 const > const m_cslSlope;

  /// A reference to the ArrayView holding the shape parameter for each element.
  arrayView1d< real64 const > const m_shapeParameter;

  /// A reference to the ArrayView holding the new preconsolidation pressure for each integration point
  arrayView2d< real64 > const m_newPreConsolidationPressure;

  /// A reference to the ArrayView holding the old preconsolidation presure for each integration point
  arrayView2d< real64 > const m_oldPreConsolidationPressure;

};


GEOS_HOST_DEVICE
inline
void DelftEggUpdates::evaluateYield( real64 const p,
                                     real64 const q,
                                     real64 const pc,
                                     real64 const M,
                                     real64 const alpha,
                                     real64 const Cc,
                                     real64 const Cr,
                                     real64 const bulkModulus,
                                     real64 const mu,
                                     real64 & f,
                                     real64 & df_dp,
                                     real64 & df_dq,
                                     real64 & df_dpc,
                                     real64 & df_dp_dve,
                                     real64 & df_dq_dse ) const
{
  real64 const c = alpha / (alpha+1.0) * pc;
  real64 a = alpha;
  real64 pa = pc;
  real64 factor = 1.0;

  if( p >= c ) // Use MCC
  {
    a = 1.0;
    factor = 2.0 * alpha / (alpha+1.0);
    pa = factor * pc;
  }
  real64 alphaTerm = 2.0 * a * a * a / (a+1.0);
  df_dp = -alphaTerm * pa + 2.0 * a * a * p;
  df_dq = 2.0 * q / (M * M);
  df_dpc = 2.0 * a * a * (a-1.0) / (a+1.0) * pc - alphaTerm * p * factor;
  real64 dpc_dve = -1.0 / (Cc-Cr) * pc;
  df_dp_dve = 2.0 * a * a * bulkModulus + alphaTerm * dpc_dve * factor;
  df_dq_dse = 2.0 / (M * M) * 3.0 * mu;

  f = q * q / (M * M)- a * a * p *(2.0 * a / (a+1.0) * pa - p) + a * a * (a-1.0) / (a+1.0) * pc * pc;

}


GEOS_HOST_DEVICE
inline
void DelftEggUpdates::smallStrainUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const & timeIncrement,
                                         real64 const ( &strainIncrement )[6],
                                         real64 ( & stress )[6],
                                         real64 ( & stiffness )[6][6] ) const
{

  // Rename variables for easier implementation
  real64 const oldPc  = m_oldPreConsolidationPressure[k][q];   //pre-consolidation pressure
  real64 const mu     = m_shearModulus[k];
  real64 const bulkModulus     = m_bulkModulus[k];

  real64 const M      = m_cslSlope[k];
  real64 const Cr     = m_recompressionIndex[k];
  real64 const Cc     = m_virginCompressionIndex[k];
  real64 const alpha  = m_shapeParameter[k];

  real64 pc = oldPc;

  // elastic predictor (assume strainIncrement is all elastic)

  ElasticIsotropicUpdates::smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );

  if( m_disableInelasticity )
  {
    return;
  }

  // check yield function F <= 0

  real64 trialP;
  real64 trialQ;
  real64 eps_v_trial;
  real64 eps_s_trial;
  real64 deviator[6];

  twoInvariant::stressDecomposition( stress,
                                     trialP,
                                     trialQ,
                                     deviator );

  real64 yield, df_dp, df_dq, df_dpc, df_dp_dve, df_dq_dse;
  evaluateYield( trialP, trialQ, pc, M, alpha, Cc, Cr, bulkModulus, mu, yield, df_dp, df_dq, df_dpc, df_dp_dve, df_dq_dse );

  if( yield < 1e-9 ) // elasticity
  {
    return;
  }

  // else, plasticity (trial stress point lies outside yield surface)

  eps_v_trial = trialP/bulkModulus;

  eps_s_trial = trialQ/3.0/mu;

  real64 solution[3] = {}, residual[3] = {}, delta[3] = {};
  real64 jacobian[3][3] = {{}}, jacobianInv[3][3] = {{}};

  solution[0] = eps_v_trial; // initial guess for elastic volumetric strain
  solution[1] = eps_s_trial; // initial guess for elastic deviatoric strain
  solution[2] = 1e-10;       // initial guess for plastic multiplier

  real64 norm, normZero = 1e30;
  integer cuts = 0;
  integer const maxCuts = 5;
  real64 normOld = normZero;

  // begin Newton loop

  for( localIndex iter=0; iter<20; ++iter )
  {
    trialP = solution[0] * bulkModulus;
    trialQ = 3.0 * mu * solution[1];
    pc = oldPc * std::exp( -1.0 / (Cc-Cr) * (eps_v_trial-solution[0]));

    evaluateYield( trialP, trialQ, pc, M, alpha, Cc, Cr, bulkModulus, mu, yield, df_dp, df_dq, df_dpc, df_dp_dve, df_dq_dse );
    real64 dpc_dve = -1.0 / (Cc-Cr) * pc;

    real64 scale = 1.0 / (mu*mu);
    // assemble residual system
    residual[0] = solution[0] - eps_v_trial + solution[2]*df_dp;   // strainElasticDev - strainElasticTrialDev + dlambda*dG/dPQ = 0
    residual[1] = solution[1] - eps_s_trial + solution[2]*df_dq;         // strainElasticVol - strainElasticTrialVol + dlambda*dG/dQ = 0
    residual[2] = yield * scale;      // F = 0


    // check for convergence

    norm = LvArray::tensorOps::l2Norm< 3 >( residual );
    if( iter==0 )
    {
      normZero = norm;
      normOld = norm;
    }

    if( norm < 1e-12*(normZero+1.0))
    {
      break;
    }
    else if( iter > 0 && norm>normOld && cuts<maxCuts ) //linesearch
    {
      cuts++;
      iter--;
      LvArray::tensorOps::scale< 3 >( delta, 0.5 );
      solution[0] += delta[0];
      solution[1] += delta[1];
      solution[2] += delta[2];
      normOld = norm;
    }
    else
    {
      // solve Newton system
      normOld = norm;
      cuts=0;
      real64 dp_dve = bulkModulus;
      real64 dq_dse = 3.0 * mu;

      jacobian[0][0] = 1.0 + solution[2] * df_dp_dve;
      jacobian[0][2] = df_dp;
      jacobian[1][1] = 1.0 + solution[2] * df_dq_dse;
      jacobian[1][2] = df_dq;
      jacobian[2][0] = (dp_dve * df_dp - dpc_dve * df_dpc) * scale;
      jacobian[2][1] = (dq_dse * df_dq) * scale;
      jacobian[2][2] = 0.0;

      LvArray::tensorOps::invert< 3 >( jacobianInv, jacobian );
      LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( delta, jacobianInv, residual );

      for( localIndex i=0; i<3; ++i )
      {
        solution[i] -= delta[i];
      }
    }
  }

  // re-construct stress = P*eye + sqrt(2/3)*Q*nhat

  twoInvariant::stressRecomposition( trialP,
                                     trialQ,
                                     deviator,
                                     stress );

  // construct consistent tangent stiffness

  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0.0 );
  real64 BB[2][2] = {{}};

  real64 const c = alpha/(alpha+1.)*pc;
  real64 dpc_dve = -1./(Cc-Cr) * pc;
  real64 df_dp_depsv;
  real64 factor;
  if( trialP >= c )   // Use MCC
  {
    factor = 2.0 * alpha / (alpha+1.0);
    df_dpc = -factor * trialP;
    df_dp_depsv = factor * dpc_dve;
  }
  else
  {
    factor = 2.0 * alpha * alpha * alpha / (alpha+1.0);
    df_dpc = 2.0 * alpha * alpha * (alpha-1.0) /(alpha+1.0) * pc - factor * trialP;
    df_dp_depsv = factor * dpc_dve;
  }

  real64 a1 = 1.0 + solution[2] * df_dp_depsv;
  real64 a2 = -df_dpc * dpc_dve;
  real64 scale = 1.0 / (mu * mu); //add scaling factor to improve convergence
  BB[0][0] = bulkModulus * (a1*jacobianInv[0][0] + a2 * jacobianInv[0][2] * scale);
  BB[0][1] = bulkModulus * jacobianInv[0][1];
  BB[1][0] = 3.0 * mu * (a1 * jacobianInv[1][0] + a2 * jacobianInv[1][2] * scale);
  BB[1][1] = 3.0 * mu * jacobianInv[1][1];

  real64 c1;

  if( eps_s_trial<1e-15 ) // confirm eps_s_trial != 0
  {
    c1 = 2.0 * mu;
  }
  else
  {
    c1 = 2.0 * trialQ/(3.0 * eps_s_trial);
  }

  real64 c2 = BB[0][0] - c1 / 3.0;
  real64 c3 = std::sqrt( 2./3. ) * BB[0][1];
  real64 c4 = std::sqrt( 2./3. ) * BB[1][0];
  real64 c5 = 2./3. * BB[1][1] - c1;

  real64 identity[6];

  for( localIndex i=0; i<6; ++i )
  {
    for( localIndex j=0; j<6; ++j )
    {
      stiffness[i][j] =  0.0;
    }
  }

  for( localIndex i=0; i<3; ++i )
  {
    stiffness[i][i] = c1;
    stiffness[i+3][i+3] = 0.5 * c1;
    identity[i] = 1.0;
    identity[i+3] = 0.0;
  }

  for( localIndex i=0; i<6; ++i )
  {
    for( localIndex j=0; j<6; ++j )
    {
      stiffness[i][j] +=   c2 * identity[i] * identity[j]
                         + c3 * identity[i] * deviator[j]
                         + c4 * deviator[i] * identity[j]
                         + c5 * deviator[i] * deviator[j];
    }
  }

  // remember history variables before returning
  m_newPreConsolidationPressure[k][q] = pc;

  // save new stress and return
  saveStress( k, q, stress );
  return;
}


GEOS_HOST_DEVICE
inline
void DelftEggUpdates::smallStrainUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const & timeIncrement,
                                         real64 const ( &strainIncrement )[6],
                                         real64 ( & stress )[6],
                                         DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness.m_c );
}



/**
 * @class DelftEgg
 *
 * DelftEgg  material model.
 */
class DelftEgg : public ElasticIsotropic
{
public:

  /// @typedef Alias for DelftEggUpdates
  using KernelWrapper = DelftEggUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  DelftEgg( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~DelftEgg() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "DelftEgg";

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
    /// string/key for default recompression index
    static constexpr char const * defaultRecompressionIndexString() { return "defaultRecompressionIndex"; }

    /// string/key for virgin compression index
    static constexpr char const * defaultVirginCompressionIndexString() { return "defaultVirginCompressionIndex"; }

    /// string/key for default slope of the critical state line
    static constexpr char const * defaultCslSlopeString() { return "defaultCslSlope"; }

    /// string/key for default shape parameter
    static constexpr char const * defaultShapeParameterString() { return "defaultShapeParameter"; }

    /// string/key for default pre-consolidation pressure
    static constexpr char const * defaultPreConsolidationPressureString() { return "defaultPreConsolidationPressure"; }

    /// string/key for recompression index
    static constexpr char const * recompressionIndexString() { return "recompressionIndex"; }

    /// string/key for virgin compression index
    static constexpr char const * virginCompressionIndexString() { return "virginCompressionIndex"; }

    /// string/key for slope of the critical state line
    static constexpr char const * cslSlopeString() { return "cslSlope"; }

    /// string/key for shape parameter of the yield surface
    static constexpr char const * shapeParameterString() { return "shapeParameter"; }

    /// string/key for new pre-consolidation pressure
    static constexpr char const * newPreConsolidationPressureString() { return "preConsolidationPressure"; }

    /// string/key for old pre-consolidation pressure
    static constexpr char const * oldPreConsolidationPressureString() { return "oldPreConsolidationPressure"; }
  };

  /**
   * @brief Create a instantiation of the DelftEggUpdate class that refers to the data in this.
   * @return An instantiation of DelftEggUpdate.
   */
  DelftEggUpdates createKernelUpdates() const
  {
    return DelftEggUpdates( m_recompressionIndex,
                            m_virginCompressionIndex,
                            m_cslSlope,
                            m_shapeParameter,
                            m_newPreConsolidationPressure,
                            m_oldPreConsolidationPressure,
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
                          m_recompressionIndex,
                          m_virginCompressionIndex,
                          m_cslSlope,
                          m_shapeParameter,
                          m_newPreConsolidationPressure,
                          m_oldPreConsolidationPressure,
                          m_bulkModulus,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_disableInelasticity );
  }


protected:
  virtual void postInputInitialization() override;

  /// Material parameter: The default value of the recompression index
  real64 m_defaultRecompressionIndex;

  /// Material parameter: The default value of the virgin compression index
  real64 m_defaultVirginCompressionIndex;

  /// Material parameter: The default value of the slope of the critical state line
  real64 m_defaultCslSlope;

  /// Material parameter: The default value of the shape parameter of the yield surface
  real64 m_defaultShapeParameter;

  /// Material parameter: The default value of the preconsolidation pressure
  real64 m_defaultPreConsolidationPressure;

  /// Material parameter: The recompression index for each element
  array1d< real64 > m_recompressionIndex;

  /// Material parameter: The virgin compression index for each element
  array1d< real64 > m_virginCompressionIndex;

  /// Material parameter: The slope of the critical state line for each element
  array1d< real64 > m_cslSlope;

  /// Material parameter: The shape parameter of the yield surface for each element
  array1d< real64 > m_shapeParameter;

  /// State variable: The current preconsolidation pressure for each quadrature point
  array2d< real64 > m_newPreConsolidationPressure;

  /// State variable: The previous preconsolidation pressure for each quadrature point
  array2d< real64 > m_oldPreConsolidationPressure;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_DELFTEGG_HPP_ */
