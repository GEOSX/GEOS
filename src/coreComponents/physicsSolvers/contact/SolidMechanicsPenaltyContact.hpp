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
 * @file SolidMechanicsPenaltyContact.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSPENALTYCONTACT_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSPENALTYCONTACT_HPP_

#include "physicsSolvers/contact/ContactSolverBase.hpp"
#include "../../linearAlgebra/DofManager.hpp"
#include "../../common/DataTypes.hpp"

namespace geos
{

class SolidMechanicsLagrangianFEM;

class SolidMechanicsPenaltyContact : public ContactSolverBase
{
public:

  SolidMechanicsPenaltyContact( const string & name,
                                Group * const parent );

  ~SolidMechanicsPenaltyContact() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SolidMechanicsPenaltyContact";
  }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "PenaltyContact"; }

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true ) override final;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  void assembleContact( DomainPartition & domain,
                        DofManager const & dofManager,
                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                        arrayView1d< real64 > const & localRhs );

protected:

private:

  struct viewKeyStruct : ContactSolverBase::viewKeyStruct
  {
    constexpr static char const * contactPenaltyStiffnessString() { return "contactPenaltyStiffness"; }
  };

  real64 m_contactPenaltyStiffness;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSPENALTYCONTACT_HPP_ */
