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
 * @file FiniteElementBase.hpp
 */

#if defined(GEOS_USE_DEVICE)
#define CALC_FEM_SHAPE_IN_KERNEL
#endif



#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"
#include "finiteElement/PDEUtilities.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/EdgeManager.hpp"
#include "mesh/FaceManager.hpp"

namespace geos
{
namespace finiteElement
{

/**
 * @brief Base class for FEM element implementations.
 */
class FiniteElementBase
{
public:

  /// Default Constructor
  FiniteElementBase() = default;


  /// Number of sampling points.
  constexpr static int numSamplingPointsPerDirection = 10;

  /**
   * @brief Copy Constructor
   * @param source The object to copy.
   */
  GEOS_HOST_DEVICE
  FiniteElementBase( FiniteElementBase const & source ):
#ifdef CALC_FEM_SHAPE_IN_KERNEL
    m_viewGradN(),
    m_viewDetJ()
#else
    m_viewGradN( source.m_viewGradN ),
    m_viewDetJ( source.m_viewDetJ )
#endif
  {
    GEOS_UNUSED_VAR( source ); // suppress warning when CALC_FEM_SHAPE_IN_KERNEL is defined
  }

  /// Default Move constructor
  FiniteElementBase( FiniteElementBase && ) = default;

  /**
   * @brief Deleted copy assignment operator
   * @return deleted
   */
  FiniteElementBase & operator=( FiniteElementBase const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   * @return deleted
   */
  FiniteElementBase & operator=( FiniteElementBase && ) = delete;

  /**
   * @brief Destructor
   */
  GEOS_HOST_DEVICE
  virtual ~FiniteElementBase()
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables allocated on the stack.
   *
   * Contains variables that will be allocated on the stack. Used only by Virtual Element classes to
   * hold the computed projections of basis functions
   */
  struct StackVariables
  {
    /**
     * Default constructor
     */
    GEOS_HOST_DEVICE
    StackVariables()
    {}
  };

  /**
   * @struct MeshData
   * @brief Variables used to initialize the class.
   */
  template< typename SUBREGION_TYPE >
  struct MeshData
  {
    /**
     * Constructor
     */
    MeshData()
    {}
  };


  /**
   * @brief Method to fill a MeshData object.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param cellSubRegion The cell sub-region for which the element has to be initialized.
   * @param meshData MeshData struct to be filled.
   */
  template< typename SUBREGION_TYPE >
  static void fillMeshData( NodeManager const & nodeManager,
                            EdgeManager const & edgeManager,
                            FaceManager const & faceManager,
                            SUBREGION_TYPE const & cellSubRegion,
                            MeshData< SUBREGION_TYPE > & meshData )
  {
    GEOS_UNUSED_VAR( nodeManager,
                     edgeManager,
                     faceManager,
                     cellSubRegion,
                     meshData );
  }

  /**
   * @brief Abstract initialization method.
   * @details It calls the fillMeshData method of the specific element implementation.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param cellSubRegion The cell sub-region for which the element has to be initialized.
   * @param meshData The struct to be filled according to the @p LEAF class needs.
   */
  template< typename LEAF, typename SUBREGION_TYPE >
  static void initialize( NodeManager const & nodeManager,
                          EdgeManager const & edgeManager,
                          FaceManager const & faceManager,
                          SUBREGION_TYPE const & cellSubRegion,
                          typename LEAF::template MeshData< SUBREGION_TYPE > & meshData
                          )
  {
    LEAF::template fillMeshData< SUBREGION_TYPE >( nodeManager,
                                                   edgeManager,
                                                   faceManager,
                                                   cellSubRegion,
                                                   meshData );
  }


  /**
   * @brief Empty setup method.
   * @param cellIndex The index of the cell with respect to the cell sub region.
   * @param meshData MeshData struct filled by @ref fillMeshData.
   * @param stack Object that holds stack variables.
   */
  template< typename SUBREGION_TYPE >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void setupStack( localIndex const & cellIndex,
                          MeshData< SUBREGION_TYPE > const & meshData,
                          StackVariables & stack )
  {
    GEOS_UNUSED_VAR( cellIndex,
                     meshData,
                     stack );

  }

  /**
   * @brief Abstract setup method, possibly computing cell-dependent properties.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param cellIndex The index of the cell with respect to the cell sub region to which the element
   * has been initialized previously (see @ref initialize).
   * @param meshData A MeshData object previously filled.
   * @param stack Object that holds stack variables.
   */
  template< typename LEAF, typename SUBREGION_TYPE >
  GEOS_HOST_DEVICE
  void setup( localIndex const & cellIndex,
              typename LEAF::template MeshData< SUBREGION_TYPE > const & meshData,
              typename LEAF::StackVariables & stack ) const
  {
    LEAF::setupStack( cellIndex, meshData, stack );
  }

  /**
   * @brief Virtual getter for the number of quadrature points per element.
   * @return The number of quadrature points per element.
   */
  GEOS_HOST_DEVICE
  virtual localIndex getNumQuadraturePoints() const = 0;

  /**
   * @brief Virtual getter for the number of support points per element.
   * @return The number of support points per element.
   */
  GEOS_HOST_DEVICE
  virtual localIndex getNumSupportPoints() const = 0;

  /**
   * @brief An helper struct to determine the function space.
   * @tparam N The number of components per support point (i.e., 1 if
   *   scalar variable, 3 if vector variable)
   */
  template< int N >
  struct FunctionSpaceHelper
  {};

  /**
   * @brief Getter for the function space.
   * @tparam The number of components per support point (i.e., 1 if
   *   scalar variable, 3 if vector variable)
   * @return The function space.
   */
  template< int N >
  GEOS_HOST_DEVICE
  constexpr static PDEUtilities::FunctionSpace getFunctionSpace();

  /**
   * @brief Getter for the number of support points per element.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param stack Stack variables created by a call to @ref setup.
   * @return The number of support points per element.
   */
  template< typename LEAF >
  GEOS_HOST_DEVICE
  localIndex numSupportPoints( typename LEAF::StackVariables const & stack ) const
  {
    return LEAF::getNumSupportPoints( stack );
  }

  /**
   * @brief Get the maximum number of support points for this element.
   * @details This should be used to know the size of pre-allocated objects whose size depend on the
   * number of support points.
   * @return The number of maximum support points for this element.
   */
  GEOS_HOST_DEVICE
  virtual localIndex getMaxSupportPoints() const = 0;

  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X Array of coordinates as the reference for the gradients.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function calls the function to calculate shape function gradients.
   */
  template< typename LEAF >
  GEOS_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   real64 const (&X)[LEAF::maxSupportPoints][3],
                   real64 ( &gradN )[LEAF::maxSupportPoints][3] ) const;

  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X Array of coordinates as the reference for the gradients.
   * @param stack Stack variables relative to the element @p k created by a call to @ref setup.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function calls the function to calculate shape function gradients.
   */
  template< typename LEAF >
  GEOS_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   real64 const (&X)[LEAF::maxSupportPoints][3],
                   typename LEAF::StackVariables const & stack,
                   real64 ( &gradN )[LEAF::maxSupportPoints][3] ) const;

  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X dummy variable.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function returns pre-calculated shape function gradients.
   */
  template< typename LEAF >
  GEOS_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   int const X,
                   real64 ( &gradN )[LEAF::maxSupportPoints][3] ) const;
  /**
   * @brief Get the shape function gradients.
   * @tparam LEAF Type of the derived finite element implementation.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param X dummy variable.
   * @param stack Stack variables relative to the element @p k created by a call to @ref setup.
   * @param gradN Return array of the shape function gradients.
   * @return The determinant of the Jacobian transformation matrix.
   *
   * This function returns pre-calculated shape function gradients.
   */
  template< typename LEAF >
  GEOS_HOST_DEVICE
  real64 getGradN( localIndex const k,
                   localIndex const q,
                   int const X,
                   typename LEAF::StackVariables const & stack,
                   real64 ( &gradN )[LEAF::maxSupportPoints][3] ) const;


  /**
   * @brief Empty method, here for compatibility with methods that require a stabilization of the
   * grad-grad bilinear form.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @tparam MAXSUPPORTPOINTS Maximum number of support points allowed for this element.
   * @tparam UPPER If true only the upper triangular part of @p matrix is modified.
   * @param stack Stack variables as filled by @ref setupStack.
   * @param matrix The matrix that needs to be stabilized.
   * @param scaleFactor Scaling of the stabilization matrix.
   */
  template< localIndex NUMDOFSPERTRIALSUPPORTPOINT, localIndex MAXSUPPORTPOINTS, bool UPPER >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void addGradGradStabilization( StackVariables const & stack,
                                        real64 ( & matrix )
                                        [MAXSUPPORTPOINTS * NUMDOFSPERTRIALSUPPORTPOINT]
                                        [MAXSUPPORTPOINTS * NUMDOFSPERTRIALSUPPORTPOINT],
                                        real64 const & scaleFactor )
  {
    GEOS_UNUSED_VAR( stack,
                     matrix,
                     scaleFactor );
  }


  /**
   * @brief Add stabilization of grad-grad bilinear form to input matrix.
   * @tparam LEAF Type of the derived finite element implementation.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @tparam UPPER If true only the upper triangular part of @p matrix is modified.
   * @param stack Stack variables created by a call to @ref setup.
   * @param matrix The input matrix to which values have to be added.
   * @param scaleFactor Optional scaling of the stabilization matrix. Defaults to 1.0.
   */
  template< typename LEAF, localIndex NUMDOFSPERTRIALSUPPORTPOINT, bool UPPER = false >
  GEOS_HOST_DEVICE
  void addGradGradStabilizationMatrix( typename LEAF::StackVariables const & stack,
                                       real64 ( & matrix )
                                       [LEAF::maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT]
                                       [LEAF::maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT],
                                       real64 const scaleFactor = 1.0 ) const
  {
    LEAF::template addGradGradStabilization< NUMDOFSPERTRIALSUPPORTPOINT,
                                             LEAF::maxSupportPoints,
                                             UPPER >( stack,
                                                      matrix,
                                                      scaleFactor );
  }

  /**
   * @brief Empty method, here for compatibility with methods that require a stabilization of the
   * grad-grad bilinear form.
   * @details This method is intended to be used with @p targetVector being the residual and @p dofs
   * being the degrees of freedom of the previous solution.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @param stack Stack variables as filled by @ref setupStack.
   * @param dofs The degrees of freedom of the function where the stabilization operator has to be
   * evaluated.
   * @param targetVector The input vector to which values have to be added, seen in chunks of length
   * @p NUMDOFSPERTRIALSUPPORTPOINT.
   * @param scaleFactor Scaling of the stabilization matrix.
   */
  template< localIndex NUMDOFSPERTRIALSUPPORTPOINT, localIndex MAXSUPPORTPOINTS >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void addEvaluatedGradGradStabilization( StackVariables const & stack,
                                                 real64 const ( &dofs )[MAXSUPPORTPOINTS][NUMDOFSPERTRIALSUPPORTPOINT],
                                                 real64 ( & targetVector )[MAXSUPPORTPOINTS][NUMDOFSPERTRIALSUPPORTPOINT],
                                                 real64 const scaleFactor )
  {
    GEOS_UNUSED_VAR( stack );
    GEOS_UNUSED_VAR( dofs );
    GEOS_UNUSED_VAR( targetVector );
    GEOS_UNUSED_VAR( scaleFactor );
  }

  /**
   * @brief Add a grad-grad stabilization operator evaluated at a provided vector of dofs to input
   * vector.
   * @details This method is used to modify a residual consistently when the jacobian includes a
   * stabilization term.
   * @tparam LEAF Type of the derived finite element implementation.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @param stack Stack variables created by a call to @ref setup.
   * @param dofs The vector of dofs to evaluate the stabilization.
   * @param targetVector The input vector to which values have to be added, seen in chunks of length
   * @p NUMDOFSPERTRIALSUPPORTPOINT.
   * @param scaleFactor Optional scaling of the stabilization matrix. Defaults to 1.0.
   */
  template< typename LEAF, localIndex NUMDOFSPERTRIALSUPPORTPOINT >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void
  addEvaluatedGradGradStabilizationVector( typename LEAF::StackVariables const & stack,
                                           real64 const ( &dofs )[LEAF::maxSupportPoints]
                                           [NUMDOFSPERTRIALSUPPORTPOINT],
                                           real64 ( & targetVector )[LEAF::maxSupportPoints]
                                           [NUMDOFSPERTRIALSUPPORTPOINT],
                                           real64 const scaleFactor = 1.0 ) const
  {
    LEAF::template
    addEvaluatedGradGradStabilization< NUMDOFSPERTRIALSUPPORTPOINT, LEAF::maxSupportPoints >( stack,
                                                                                              dofs,
                                                                                              targetVector,
                                                                                              scaleFactor );
  }

  /**
   * @name Value Operator Functions
   */
  ///@{

  /**
   * @brief Compute the interpolated value of a variable.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @param N Array (for each support point) of shape function values at the
   *   coordinate the variable is to be interpolated.
   * @param var Array of variable values for each support point.
   * @param value The interpolated value of @p var.
   *
   * This is the standard finite element interpolation operator of a discrete
   * variable defined at the support points.
   * The operator is expressed as:
   * \f[
   * value  = \sum_a^{numSupport} \left ( N_a var_a \right ),
   * \f]

   * @note The shape function values @p N must be evaluated prior to calling this
   * function.
   *
   */
  template< int NUM_SUPPORT_POINTS >
  GEOS_HOST_DEVICE
  static
  void value( real64 const (&N)[NUM_SUPPORT_POINTS],
              real64 const (&var)[NUM_SUPPORT_POINTS],
              real64 & value );

  /**
   * @brief Compute the interpolated value of a vector variable.
   * @tparam NUM_COMPONENTS Number of components for the vector variable.
   * @copydoc value
   */
  template< int NUM_SUPPORT_POINTS,
            int NUM_COMPONENTS >
  GEOS_HOST_DEVICE
  static
  void value( real64 const (&N)[NUM_SUPPORT_POINTS],
              real64 const (&var)[NUM_SUPPORT_POINTS][NUM_COMPONENTS],
              real64 ( &value )[NUM_COMPONENTS] );

  ///@}

  /**
   * @name Gradient Operator Functions
   */
  ///@{

  /**
   * @brief Calculate the symmetric gradient of a vector valued support field
   *   at a point using the stored basis function gradients for all support
   *   points.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   * @param gradN The basis function gradients at a point in the element.
   * @param var The vector valued support field that the gradient operator will
   *  be applied to.
   * @param gradVar The symmetric gradient in Voigt notation.
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad^s_{ij}  = \frac{1}{2} \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ai} + \frac{\partial N_a}{\partial X_i}
   * var_{aj}\right ),
   * \f]
   *
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOS_HOST_DEVICE
  static void symmetricGradient( GRADIENT_TYPE const & gradN,
                                 real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                 real64 ( &gradVar )[6] );

  /**
   * @brief Calculate the trace of the symmetric gradient of a vector valued support
   *   field (i.e. the volumetric strain for the displacement field) at a point using
   *   the stored basis function gradients for all support points.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   * @param gradN The basis function gradients at a point in the element.
   * @param var The vector valued support field that the gradient operator will
   *  be applied to.
   * @return The trace of the symetric gradient tensor.
   *
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOS_HOST_DEVICE
  static real64 symmetricGradientTrace( GRADIENT_TYPE const & gradN,
                                        real64 const (&var)[NUM_SUPPORT_POINTS][3] );

  /**
   * @brief Calculate the gradient of a scalar valued support field at a point
   *   using the input basis function gradients.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   *   function gradients.
   * @param gradN The basis function gradients at a point in the element.
   * @param var The vector valued support field that the gradient operator will
   *  be applied to.
   * @param gradVar The  gradient.
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad_{j}  = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{a}\right ),
   * \f]
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOS_HOST_DEVICE
  static void gradient( GRADIENT_TYPE const & gradN,
                        real64 const (&var)[NUM_SUPPORT_POINTS],
                        real64 ( &gradVar )[3] );

  /**
   * @brief Calculate the gradient of a vector valued support field at a point
   *   using the input basis function gradients.
   * @copydoc gradient
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad_{ij}  = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ai}\right ),
   * \f]
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOS_HOST_DEVICE
  static void gradient( GRADIENT_TYPE const & gradN,
                        real64 const (&var)[NUM_SUPPORT_POINTS][3],
                        real64 ( &gradVar )[3][3] );
  ///@}

  /**
   * @name Multi-Operator Functions
   */
  ///@{

  /**
   * @brief Calculate the value and gradient of a scalar valued support field
   *   at a point using the input basis function gradients.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   * @param N Array (for each support point) of shape function values at the
   *   coordinate the variable is to be interpolated.
   * @param gradN The basis function gradients at a point in the element.
   * @param var The vector valued support field that the gradient operator will
   *  be applied to.
   * @param value The value at the point for which N was specified.
   * @param gradVar The gradient at the point for which gradN was specified.
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOS_HOST_DEVICE
  static void valueAndGradient( real64 const (&N)[NUM_SUPPORT_POINTS],
                                GRADIENT_TYPE const & gradN,
                                real64 const (&var)[NUM_SUPPORT_POINTS],
                                real64 & value,
                                real64 ( &gradVar )[3] );

  ///@}

  /**
   * @name Scattering Operator Functions
   *
   * These functions take quadrature data and map it to the support points
   * through some operator.
   */
  ///@{

  /**
   * @brief Inner product of each basis function gradient with a rank-2
   *   symmetric tensor.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   *   function gradients.
   * @param gradN The basis function gradients at a point in the element.
   * @param var_detJxW The rank-2 tensor at @p q scaled by J*W.
   * @param R The vector at each support point which will hold the result from
   *   the tensor contraction.
   *
   * More precisely, the operator is defined as:
   *
   * \f[
   * R_i = \sum_a^{nSupport} \left( \frac{\partial N_a}{\partial X_j} var_{ij}\right),
   * \f]
   *
   * where \f$\frac{\partial N_a}{\partial X_j}\f$ is the basis function gradient,
   *   \f$var_{ij}\f$ is the rank-2 symmetric tensor.
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOS_HOST_DEVICE
  static void plusGradNajAij( GRADIENT_TYPE const & gradN,
                              real64 const (&var_detJxW)[6],
                              real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  /**
   * @copydoc plusGradNajAij
   * @brief Inner product of each basis function gradient with a rank-2
   *   tensor.
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOS_HOST_DEVICE
  static void plusGradNajAij( GRADIENT_TYPE const & gradN,
                              real64 const (&var_detJxW)[3][3],
                              real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  /**
   * @brief Product of each shape function with a vector forcing term.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @param N The shape function value at a predetermined coordinate in the element.
   * @param forcingTerm_detJxW A vector scaled by detJxW
   * @param R The vector at each support point which will hold the result from
   *   the tensor contraction.
   */
  template< int NUM_SUPPORT_POINTS >
  GEOS_HOST_DEVICE
  static void plusNaFi( real64 const (&N)[NUM_SUPPORT_POINTS],
                        real64 const (&forcingTerm_detJxW)[3],
                        real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  /**
   * @brief Inner product of each basis function gradient with a rank-2
   *   symmetric tensor added to the product each shape function with a vector.
   * @tparam NUM_SUPPORT_POINTS The number of support points for the element.
   * @tparam GRADIENT_TYPE The type of the array object holding the shape
   *   function gradients.
   * @param gradN The basis function gradients at a point in the element.
   * @param var_detJxW The rank-2 symmetric tensor at @p q scaled by J*W.
   * @param N The shape function value at a predetermined coordinate in the element.
   * @param forcingTerm_detJxW A vector scaled by detJxW
   * @param R The vector at each support point which will hold the result from
   *   the tensor contraction.
   *
   * \f[
   * R_i = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ij} + N_a f_i \right ),
   * \f]
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOS_HOST_DEVICE
  static void plusGradNajAijPlusNaFi( GRADIENT_TYPE const & gradN,
                                      real64 const (&var_detJxW)[3][3],
                                      real64 const (&N)[NUM_SUPPORT_POINTS],
                                      real64 const (&forcingTerm_detJxW)[3],
                                      real64 ( &R )[NUM_SUPPORT_POINTS][3] );

  /**
   * @brief Inner product of each basis function gradient with a rank-2
   *   tensor added to the product each shape function with a vector.
   * @copydoc plusGradNajAijPlusNaFi
   */
  template< int NUM_SUPPORT_POINTS,
            typename GRADIENT_TYPE >
  GEOS_HOST_DEVICE
  static void plusGradNajAijPlusNaFi( GRADIENT_TYPE const & gradN,
                                      real64 const (&var_detJxW)[6],
                                      real64 const (&N)[NUM_SUPPORT_POINTS],
                                      real64 const (&forcingTerm_detJxW)[3],
                                      real64 ( &R )[NUM_SUPPORT_POINTS][3] );


  /**
   * @brief Sets m_viewGradN equal to an input view.
   * @param source The view to assign to m_viewGradN.
   */
  void setGradNView( arrayView4d< real64 const > const & source )
  {
    GEOS_ERROR_IF_NE_MSG( source.size( 1 ),
                          getNumQuadraturePoints(),
                          "2nd-dimension of gradN array does not match number of quadrature points" );
    GEOS_ERROR_IF_NE_MSG( source.size( 2 ),
                          getMaxSupportPoints(),
                          "3rd-dimension of gradN array does not match number of support points" );
    GEOS_ERROR_IF_NE_MSG( source.size( 3 ),
                          3,
                          "4th-dimension of gradN array does not match 3" );

    m_viewGradN = source;
  }

  /**
   * @brief Sets m_viewDetJ equal to an input view.
   * @param source The view to assign to m_viewDetJ.
   */
  void setDetJView( arrayView2d< real64 const > const & source )
  {
    GEOS_ERROR_IF_NE_MSG( source.size( 1 ),
                          getNumQuadraturePoints(),
                          "2nd-dimension of gradN array does not match number of quadrature points" );
    m_viewDetJ = source;
  }

  /**
   * @brief Getter for m_viewGradN
   * @return A new arrayView copy of m_viewGradN.
   */
  arrayView4d< real64 const > getGradNView() const
  {
    return m_viewGradN;
  }

  /**
   * @brief Getter for m_viewDetJ
   * @return A new arrayView copy of m_viewDetJ.
   */
  arrayView2d< real64 const > getDetJView() const
  {
    return m_viewDetJ;
  }


protected:
  /// View to potentially hold pre-calculated shape function gradients.
  arrayView4d< real64 const > m_viewGradN;

  /// View to potentially hold pre-calculated weighted jacobian transformation
  /// determinants.
  arrayView2d< real64 const > m_viewDetJ;
};

/// @cond Doxygen_Suppress

//*************************************************************************************************
//***** Definitions *******************************************************************************
//*************************************************************************************************

template<>
struct FiniteElementBase::FunctionSpaceHelper< 1 >
{
  GEOS_HOST_DEVICE
  constexpr static PDEUtilities::FunctionSpace getFunctionSpace()
  {
    return PDEUtilities::FunctionSpace::H1;
  }
};

template<>
struct FiniteElementBase::FunctionSpaceHelper< 3 >
{
  GEOS_HOST_DEVICE
  constexpr static PDEUtilities::FunctionSpace getFunctionSpace()
  {
    return PDEUtilities::FunctionSpace::H1vector;
  }
};

template< int N >
GEOS_HOST_DEVICE
constexpr PDEUtilities::FunctionSpace FiniteElementBase::getFunctionSpace()
{
  return FunctionSpaceHelper< N >::getFunctionSpace();
}

template< typename LEAF >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    real64 const (&X)[LEAF::maxSupportPoints][3],
                                    real64 (& gradN)[LEAF::maxSupportPoints][3] ) const
{
  GEOS_UNUSED_VAR( k );
  return LEAF::calcGradN( q, X, gradN );
}

template< typename LEAF >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    real64 const (&X)[LEAF::maxSupportPoints][3],
                                    typename LEAF::StackVariables const & stack,
                                    real64 ( & gradN )[LEAF::maxSupportPoints][3] ) const
{
  GEOS_UNUSED_VAR( k );
  return LEAF::calcGradN( q, X, stack, gradN );
}

template< typename LEAF >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    int const X,
                                    real64 (& gradN)[LEAF::maxSupportPoints][3] ) const
{
  GEOS_UNUSED_VAR( X );

  LvArray::tensorOps::copy< LEAF::maxSupportPoints, 3 >( gradN, m_viewGradN[ k ][ q ] );

  return m_viewDetJ( k, q );
}

template< typename LEAF >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FiniteElementBase::getGradN( localIndex const k,
                                    localIndex const q,
                                    int const X,
                                    typename LEAF::StackVariables const & stack,
                                    real64 (& gradN)[LEAF::maxSupportPoints][3] ) const
{
  GEOS_UNUSED_VAR( X );
  GEOS_UNUSED_VAR( stack );

  LvArray::tensorOps::copy< LEAF::maxSupportPoints, 3 >( gradN, m_viewGradN[ k ][ q ] );

  return m_viewDetJ( k, q );
}

//*************************************************************************************************
//***** Interpolated Value Functions **************************************************************
//*************************************************************************************************

template< int NUM_SUPPORT_POINTS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::value( real64 const (&N)[NUM_SUPPORT_POINTS],
                               real64 const (&var)[NUM_SUPPORT_POINTS],
                               real64 & value )
{
  value = LvArray::tensorOps::AiBi< NUM_SUPPORT_POINTS >( N, var );
}

template< int NUM_SUPPORT_POINTS,
          int NUM_COMPONENTS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::value( real64 const (&N)[NUM_SUPPORT_POINTS],
                               real64 const (&var)[NUM_SUPPORT_POINTS][NUM_COMPONENTS],
                               real64 (& value)[NUM_COMPONENTS] )
{

  LvArray::tensorOps::Ri_eq_AjiBj< 3, NUM_SUPPORT_POINTS >( value, var, N );
}


//*************************************************************************************************
//***** Variable Gradient Functions ***************************************************************
//*************************************************************************************************

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::symmetricGradient( GRADIENT_TYPE const & gradN,
                                           real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                           real64 (& gradVar)[6] )
{
  gradVar[0] = gradN[0][0] * var[0][0];
  gradVar[1] = gradN[0][1] * var[0][1];
  gradVar[2] = gradN[0][2] * var[0][2];
  gradVar[3] = gradN[0][2] * var[0][1] + gradN[0][1] * var[0][2];
  gradVar[4] = gradN[0][2] * var[0][0] + gradN[0][0] * var[0][2];
  gradVar[5] = gradN[0][1] * var[0][0] + gradN[0][0] * var[0][1];

  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
  {
    gradVar[0] = gradVar[0] + gradN[a][0] * var[ a ][0];
    gradVar[1] = gradVar[1] + gradN[a][1] * var[ a ][1];
    gradVar[2] = gradVar[2] + gradN[a][2] * var[ a ][2];
    gradVar[3] = gradVar[3] + gradN[a][2] * var[ a ][1] + gradN[a][1] * var[ a ][2];
    gradVar[4] = gradVar[4] + gradN[a][2] * var[ a ][0] + gradN[a][0] * var[ a ][2];
    gradVar[5] = gradVar[5] + gradN[a][1] * var[ a ][0] + gradN[a][0] * var[ a ][1];
  }
}

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FiniteElementBase::symmetricGradientTrace( GRADIENT_TYPE const & gradN,
                                                  real64 const (&var)[NUM_SUPPORT_POINTS][3] )
{
  real64 result = gradN[0][0] * var[0][0] + gradN[0][1] * var[0][1] + gradN[0][2] * var[0][2];

  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
  {
    result = result + gradN[a][0] * var[a][0] + gradN[a][1] * var[a][1] + gradN[a][2] * var[a][2];
  }
  return result;
}

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::gradient( GRADIENT_TYPE const & gradN,
                                  real64 const (&var)[NUM_SUPPORT_POINTS],
                                  real64 (& gradVar)[3] )
{
  LvArray::tensorOps::Ri_eq_AjiBj< 3, NUM_SUPPORT_POINTS >( gradVar, gradN, var );
}

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::gradient( GRADIENT_TYPE const & gradN,
                                  real64 const (&var)[NUM_SUPPORT_POINTS][3],
                                  real64 (& gradVar)[3][3] )
{
  LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, NUM_SUPPORT_POINTS >( gradVar, var, gradN );
}



template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::valueAndGradient( real64 const (&N)[NUM_SUPPORT_POINTS],
                                          GRADIENT_TYPE const & gradN,
                                          real64 const (&var)[NUM_SUPPORT_POINTS],
                                          real64 & value,
                                          real64 (& gradVar)[3] )
{
  value = N[0] * var[0];
  for( int i = 0; i < 3; ++i )
  {
    gradVar[i] = var[0] * gradN[0][i];
  }

  for( int a=1; a<NUM_SUPPORT_POINTS; ++a )
  {
    value = value + N[a] * var[a];
    for( int i = 0; i < 3; ++i )
    {
      gradVar[i] = gradVar[i] + var[ a ] * gradN[a][i];
    }
  }
}



template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::plusGradNajAij( GRADIENT_TYPE const & gradN,
                                        real64 const (&var_detJxW)[6],
                                        real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] + var_detJxW[0] * gradN[a][0] + var_detJxW[5] * gradN[a][1] + var_detJxW[4] * gradN[a][2];
    R[a][1] = R[a][1] + var_detJxW[5] * gradN[a][0] + var_detJxW[1] * gradN[a][1] + var_detJxW[3] * gradN[a][2];
    R[a][2] = R[a][2] + var_detJxW[4] * gradN[a][0] + var_detJxW[3] * gradN[a][1] + var_detJxW[2] * gradN[a][2];
  }
}


template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::plusGradNajAij( GRADIENT_TYPE const & gradN,
                                        real64 const (&var_detJxW)[3][3],
                                        real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( R[a], var_detJxW, gradN[a] );
  }
}

template< int NUM_SUPPORT_POINTS >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::plusNaFi( real64 const (&N)[NUM_SUPPORT_POINTS],
                                  real64 const (&var_detJxW)[3],
                                  real64 ( & R )[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    LvArray::tensorOps::scaledAdd< 3 >( R[a], var_detJxW, N[a] );
  }
}


template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::plusGradNajAijPlusNaFi( GRADIENT_TYPE const & gradN,
                                                real64 const (&var_detJxW)[6],
                                                real64 const (&N)[NUM_SUPPORT_POINTS],
                                                real64 const (&forcingTerm_detJxW)[3],
                                                real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] + var_detJxW[0] * gradN[a][0] + var_detJxW[5] * gradN[a][1] + var_detJxW[4] * gradN[a][2] + forcingTerm_detJxW[0] * N[a];
    R[a][1] = R[a][1] + var_detJxW[5] * gradN[a][0] + var_detJxW[1] * gradN[a][1] + var_detJxW[3] * gradN[a][2] + forcingTerm_detJxW[1] * N[a];
    R[a][2] = R[a][2] + var_detJxW[4] * gradN[a][0] + var_detJxW[3] * gradN[a][1] + var_detJxW[2] * gradN[a][2] + forcingTerm_detJxW[2] * N[a];
  }
}

template< int NUM_SUPPORT_POINTS,
          typename GRADIENT_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FiniteElementBase::plusGradNajAijPlusNaFi( GRADIENT_TYPE const & gradN,
                                                real64 const (&var_detJxW)[3][3],
                                                real64 const (&N)[NUM_SUPPORT_POINTS],
                                                real64 const (&forcingTerm_detJxW)[3],
                                                real64 (& R)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    R[a][0] = R[a][0] + var_detJxW[0][0] * gradN[a][0] + var_detJxW[0][1] * gradN[a][1] + var_detJxW[0][2] * gradN[a][2] + forcingTerm_detJxW[0] * N[a];
    R[a][1] = R[a][1] + var_detJxW[1][0] * gradN[a][0] + var_detJxW[1][1] * gradN[a][1] + var_detJxW[1][2] * gradN[a][2] + forcingTerm_detJxW[1] * N[a];
    R[a][2] = R[a][2] + var_detJxW[2][0] * gradN[a][0] + var_detJxW[2][1] * gradN[a][1] + var_detJxW[2][2] * gradN[a][2] + forcingTerm_detJxW[2] * N[a];
  }
}
/// @endcond

}
}


/// Macro to simplify name resolution in derived classes.
#define USING_FINITEELEMENTBASE                       \
  using FiniteElementBase::value;                     \
  using FiniteElementBase::symmetricGradient;         \
  using FiniteElementBase::gradient;                  \
  using FiniteElementBase::valueAndGradient;          \
  using FiniteElementBase::plusGradNajAij;           \
  using FiniteElementBase::plusNaFi;                 \
  using FiniteElementBase::plusGradNajAijPlusNaFi;

#endif //GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_FINITEELEMENTBASE_HPP_
