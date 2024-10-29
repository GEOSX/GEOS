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
 * @file TractionBoundaryCondition.hpp
 */


#ifndef GEOS_FIELDSPECIFICATION_TRACTIONBOUNDARYCONDITION_HPP
#define GEOS_FIELDSPECIFICATION_TRACTIONBOUNDARYCONDITION_HPP

#include "FieldSpecificationBase.hpp"
#include "mesh/FaceManager.hpp"

namespace geos
{

class TableFunction;

/**
 * @class TractionBoundaryCondition
 * Holds data and methods to apply a traction boundary condition
 */
class TractionBoundaryCondition : public FieldSpecificationBase
{
public:
  /// @copydoc FieldSpecificationBase(string const &, Group *)
  TractionBoundaryCondition( string const & name, Group * parent );

  /// deleted default constructor
  TractionBoundaryCondition() = delete;

  /// default destructor
  virtual ~TractionBoundaryCondition() = default;

  /// deleted copy constructor
  TractionBoundaryCondition( TractionBoundaryCondition const & ) = delete;

  /// defaulted move constructor
  TractionBoundaryCondition( TractionBoundaryCondition && ) = default;

  /// deleted copy assignment operator
  TractionBoundaryCondition & operator=( TractionBoundaryCondition const & ) = delete;

  /// deleted move assignment operator
  TractionBoundaryCondition & operator=( TractionBoundaryCondition && ) = delete;


  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "Traction"; }


  /**
   * @brief Setup and Launche of the traction BC kernel.
   * @param time The time that should be used to evaluate any time tables.
   * @param blockLocalDofNumber Array of block local DOF numbers for the displacement.
   * @param dofRankOffset The rank offset for the DOF.
   * @param faceManager Reference to the face manager (Tractions are applied on faces)
   * @param targetSet The set of faces to apply the BC to.
   * @param localRhs The RHS of the system to add contributions to.
   */
  void launch( real64 const time,
               arrayView1d< globalIndex const > const blockLocalDofNumber,
               globalIndex const dofRankOffset,
               FaceManager const & faceManager,
               SortedArrayView< localIndex const > const & targetSet,
               arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief View keys
   */
  struct viewKeyStruct : public FieldSpecificationBase::viewKeyStruct
  {
    /// @return The key for tractionType
    constexpr static char const * tractionTypeString() { return "tractionType"; }

    /// @return The key for inputStress.
    constexpr static char const * inputStressString() { return "inputStress"; }

//    /// @return The key for the function describing the components of stress.
//    constexpr static char const * stressFunctionString() { return "stressFunctions"; }

  };

  /**
   * @brief Type of traction boundary condition.
   */
  enum class TractionType : integer
  {
    vector, ///< traction is applied to the faces as specified from the scale and direction
    normal, ///< traction is applied to the faces as a pressure specified from the product of scale and the outward face normal
    stress  ///< traction is applied to the faces as specified by the inner product of input stress and face normal
  };

protected:

  virtual void postInputInitialization() override final;

  virtual void initializePreSubGroups() override final;

  /// The type of traction to be applied, i.e. how to generate the traction.
  TractionType m_tractionType;

  /// single specified value for stress used to generate the traction if m_tractionType == stress.
  R2SymTensor m_inputStress;

//  /// names of the functions used to specify stress for the generation of tractions.
//  array1d<string> m_stressFunctionNames;
//
//  bool m_useStressFunctions;
//
//  TableFunction const * m_stressFunctions[6];

};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( TractionBoundaryCondition::TractionType,
              "vector",
              "normal",
              "stress" );

} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_TRACTIONBOUNDARYCONDITION_HPP */
