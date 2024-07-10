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

/*
 * SolidMechanicsAugmentedLagrangianContact.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSAUGMENTEDLAGRANGIANCONTACT_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSAUGMENTEDLAGRANGIANCONTACT_HPP_

#include "physicsSolvers/contact/ContactSolverBase.hpp"

namespace geos
{

class SolidMechanicsAugmentedLagrangianContact : public ContactSolverBase
{
public:
  SolidMechanicsAugmentedLagrangianContact( const string & name,
                                            Group * const parent );

  ~SolidMechanicsAugmentedLagrangianContact() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SolidMechanicsAugmentedLagrangianContact";
  }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override final;

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override final;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override final;

  virtual void implicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain ) override final;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual real64 calculateResidualNorm( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localRhs ) override;

  virtual void applySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    real64 const dt,
                                    DomainPartition & domain ) override;

  void updateState( DomainPartition & domain ) override final;

  virtual bool updateConfiguration( DomainPartition & domain ) override final;


  /**
   * @brief Loop over the finite element type on the fracture subregions of meshName and apply callback.
   * @tparam LAMBDA The callback function type
   * @param meshName The mesh name.
   * @param lambda The callback function. Take the finite element type name and
   * the list of face element of the same type.
   */
  template< typename LAMBDA >
  void forFiniteElementOnFractureSubRegions( string const & meshName, LAMBDA && lambda ) const
  {

    std::map< string,
              array1d< localIndex > > const & faceTypesToFaceElements = m_faceTypesToFaceElements.at( meshName );

    for( const auto & [finiteElementName, faceElementList] : faceTypesToFaceElements )
    {
      arrayView1d< localIndex const > const faceElemList = faceElementList.toViewConst();
      lambda( finiteElementName, faceElemList );
    }

  }

  /**
   * @brief Create the list of finite elements of the same type
   *   for each FaceElementSubRegion (Triangle or Quadrilateral).
   * @param domain The physical domain object
   */
  void createFaceTypeList( DomainPartition const & domain );

  /**
   * @brief Create the list of elements belonging to CellElementSubRegion
   *  that are enriched with the bubble basis functions
   * @param domain The physical domain object
   */
  void createBubbleCellList( DomainPartition & domain ) const;

private:

  /**
   * @brief add the number of non-zero elements induced by the coupling between
   *   nodal and bubble displacement.
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLengths the array containing the number of non-zero elements for each row
   */
  void addCouplingNumNonzeros( DomainPartition & domain,
                               DofManager & dofManager,
                               arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void addCouplingSparsityPattern( DomainPartition const & domain,
                                   DofManager const & dofManager,
                                   SparsityPatternView< globalIndex > const & pattern ) const;

  void computeTolerances( DomainPartition & domain ) const;

  /// Finite element type to face element index map
  std::map< string, std::map< string, array1d< localIndex > > > m_faceTypesToFaceElements;

  /// Finite element type to finite element object map
  std::map< string, std::unique_ptr< geos::finiteElement::FiniteElementBase > > m_faceTypeToFiniteElements;

  struct viewKeyStruct : ContactSolverBase::viewKeyStruct
  {

    constexpr static char const * normalDisplacementToleranceString() { return "normalDisplacementTolerance"; }

    constexpr static char const * normalTractionToleranceString() { return "normalTractionTolerance"; }

    constexpr static char const * slidingToleranceString() { return "slidingTolerance"; }
  };

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSAUGMENTEDLAGRANGIANCONTACT_HPP_ */
