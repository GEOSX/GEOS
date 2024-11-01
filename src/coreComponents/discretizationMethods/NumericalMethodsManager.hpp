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
 * @file NumericalMethodsManager.hpp
 */

#ifndef GEOS_DISCRETIZATIONMETHODS_NUMERICALMETHODSMANAGER_HPP_
#define GEOS_DISCRETIZATIONMETHODS_NUMERICALMETHODSMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"

namespace geos
{
namespace dataRepository
{
namespace keys
{}
}

/**
 * @class NumericalMethodsManager
 *
 * Contains the discretization methods and their components for application to the mesh.
 */
class NumericalMethodsManager : public dataRepository::Group
{
public:
  /// Deleted default constructor.
  NumericalMethodsManager() = delete;

  /**
   * @brief Constructor
   *
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  NumericalMethodsManager( string const & name, Group * const parent );

  virtual ~NumericalMethodsManager() override;

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// Contains the keys for the object names in the data repository.
  struct groupKeysStruct
  {
    /// @return Name for the basis function manager.
    static constexpr char const * basisFunctionsString() { return "BasisFunctions"; }

    /// @return Name for the quadrature rule manager.
    static constexpr char const * quadratureRulesString() { return "QuadratureRules"; }

    /// @return Name for the finite element discretization manager.
    static constexpr char const * finiteElementDiscretizationsString() { return "FiniteElements"; }

    /// @return Name for the finite volume manager.
    static constexpr char const * finiteVolumeManagerString() { return "FiniteVolume"; }
  };

  /**
   * @brief @return Returns reference to const FiniteElementDiscretizationManager m_finiteElementDiscretizationManager.
   */
  FiniteElementDiscretizationManager const & getFiniteElementDiscretizationManager() const { return m_finiteElementDiscretizationManager; }

  /**
   * @brief @return Returns reference to FiniteElementDiscretizationManager m_finiteElementDiscretizationManager.
   */
  FiniteElementDiscretizationManager & getFiniteElementDiscretizationManager()       { return m_finiteElementDiscretizationManager; }

  /**
   * @brief @return Returns reference to FiniteVolumeManager m_finiteVolumeManager.
   */
  FiniteVolumeManager & getFiniteVolumeManager()       { return m_finiteVolumeManager; }

  /**
   * @brief @return Returns reference to const FiniteVolumeManager m_finiteVolumeManager.
   */
  FiniteVolumeManager const & getFiniteVolumeManager() const { return m_finiteVolumeManager; }

private:

  /// Contains the finite element discretizations
  FiniteElementDiscretizationManager m_finiteElementDiscretizationManager;

  /// Contains the finite volume discretizations.
  FiniteVolumeManager m_finiteVolumeManager;

};

} /* namespace geos */

#endif /* GEOS_DISCRETIZATIONMETHODS_NUMERICALMETHODSMANAGER_HPP_ */
