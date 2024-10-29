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
 * @file AcousticElasticWaveEquationSEM.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEM_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEM_HPP_

#include "physicsSolvers/wavePropagation/sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEM.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "AcoustoElasticFields.hpp"
#include <tuple>

namespace geos
{

template< typename ... SOLVERS >
class CoupledWaveSolver : public SolverBase
{

public:

  /**
   * @brief main constructor for CoupledWaveSolver Objects
   * @param name the name of this instantiation of CoupledWaveSolver in the repository
   * @param parent the parent group of this instantiation of CoupledWaveSolver
   */
  CoupledWaveSolver( const string & name,
                     Group * const parent )
    : SolverBase( name, parent )
  {
    forEachArgInTuple( m_solvers, [&]( auto solver, auto idx )
    {
      using SolverType = TYPEOFPTR( solver );
      string const key = SolverType::coupledSolverAttributePrefix() + "SolverName";
      registerWrapper( key, &m_names[idx()] ).
        setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
        setInputFlag( dataRepository::InputFlags::REQUIRED ).
        setDescription( "Name of the " + SolverType::coupledSolverAttributePrefix() + " solver used by the coupled solver" );
    } );
  }

  /// deleted copy constructor
  CoupledWaveSolver( CoupledWaveSolver const & ) = delete;

  /// default move constructor
  CoupledWaveSolver( CoupledWaveSolver && ) = default;

  /// deleted assignment operator
  CoupledWaveSolver & operator=( CoupledWaveSolver const & ) = delete;

  /// deleted move operator
  CoupledWaveSolver & operator=( CoupledWaveSolver && ) = delete;

  virtual void
  postInputInitialization() override final
  {
    SolverBase::postInputInitialization();

    forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
    {
      using SolverPtr = TYPEOFREF( solver );
      using SolverType = TYPEOFPTR( SolverPtr {} );
      solver = getParent().template getGroupPointer< SolverType >( m_names[idx()] );
      GEOS_THROW_IF( solver == nullptr,
                     GEOS_FMT( "Could not find solver '{}' of type {}",
                               m_names[idx()], LvArray::system::demangleType< SolverType >() ),
                     InputError );
    } );
  }

protected:

  /// Pointers of the single-physics solvers
  std::tuple< SOLVERS *... > m_solvers;

  /// Names of the single-physics solvers
  std::array< string, sizeof...( SOLVERS ) > m_names;
};


class AcousticElasticWaveEquationSEM : public CoupledWaveSolver< AcousticWaveEquationSEM, ElasticWaveEquationSEM >
{
public:
  using Base = CoupledWaveSolver< AcousticWaveEquationSEM, ElasticWaveEquationSEM >;
  using Base::m_solvers;
  using wsCoordType = AcousticWaveEquationSEM::wsCoordType;

  enum class SolverType : integer
  {
    AcousticWaveEquationSEM = 0,
    ElasticWaveEquationSEM = 1
  };

  /// String used to form the solverName used to register solvers in CoupledWaveSolver
  static string coupledSolverAttributePrefix() { return "acousticelastic"; }

  using EXEC_POLICY = parallelDevicePolicy<  >;
  using ATOMIC_POLICY = AtomicPolicy< EXEC_POLICY >;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  /**
   * @brief main constructor for AcousticElasticWaveEquationSEM objects
   * @param name the name of this instantiation of AcousticElasticWaveEquationSEM in the repository
   * @param parent the parent group of this instantiation of AcousticElasticWaveEquationSEM
   */
  AcousticElasticWaveEquationSEM( const string & name,
                                  Group * const parent )
    : Base( name, parent )
  { }

  /// Destructor for the class
  ~AcousticElasticWaveEquationSEM() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new AcousticElasticWaveEquationSEM object through the object catalog.
   */
  static string catalogName() { return "AcousticElasticSEM"; }

  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @brief accessor for the pointer to the acoustic solver
   * @return a pointer to the acoustic solver
   */
  AcousticWaveEquationSEM * acousticSolver() const
  {
    return std::get< toUnderlying( SolverType::AcousticWaveEquationSEM ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the elastic solver
   * @return a pointer to the elastic solver
   */
  ElasticWaveEquationSEM * elasticSolver() const
  {
    return std::get< toUnderlying( SolverType::ElasticWaveEquationSEM ) >( m_solvers );
  }

  // (requires not to be private because it is called from GEOS_HOST_DEVICE method)
  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  virtual void
  cleanup( real64 const time_n,
           integer const cycleNumber,
           integer const eventCounter,
           real64 const eventProgress,
           DomainPartition & domain ) override;

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

  SortedArray< localIndex > m_interfaceNodesSet;
  arrayView1d< string const > m_acousRegions;
  arrayView1d< string const > m_elasRegions;
};

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEM_HPP_ */
