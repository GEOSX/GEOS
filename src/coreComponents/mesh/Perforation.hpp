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

/*
 * @file Perforation.hpp
 */

#ifndef GEOS_MESH_PERFORATION_HPP
#define GEOS_MESH_PERFORATION_HPP

#include "dataRepository/Group.hpp"

namespace geos
{

/**
 * @class Perforation
 *
 * This class describes a perforation with its location, well transmissibility  and corresponding well element
 */
class Perforation : public dataRepository::Group
{
public:
  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for Perforation Objects.
   * @param[in] name the name of this instantiation of Perforation in the repository
   * @param[in] parent the parent group of this instantiation of Perforation
   */
  explicit Perforation( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~Perforation() override;

  /**
   * @brief Deleted default constructor.
   */
  Perforation() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  Perforation( Perforation const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  Perforation( Perforation && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a perforation object
   */
  Perforation & operator=( Perforation const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a perforation object
   */
  Perforation & operator=( Perforation && ) = delete;

  ///@}

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the linear distance between the well head and the perforation.
   * @return the distance between the well head and the perforation
   */
  real64 const & getDistanceFromWellHead() const { return m_distanceFromHead; }

  /**
   * @brief Get the well Peaceman index at the perforation.
   * @return the well transmissibility
   */
  real64 getWellTransmissibility() const { return m_wellTransmissibility; }

  /**
   * @brief Get the well skin factor at the perforation.
   * @return the well skin factor
   */
  real64 getWellSkinFactor() const { return m_wellSkinFactor; }

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// @return String key for the linear distance from well head
    static constexpr char const *distanceFromHeadString() { return "distanceFromHead"; }
    /// @return String key for the well transmissibility at this perforation
    static constexpr char const *wellTransmissibilityString() { return "transmissibility"; }
    /// @return String key for the well skin factor at this perforation
    static constexpr char const *wellSkinFactorString() { return "skinFactor"; }
    /// ViewKey for the linear distance from well head
    dataRepository::ViewKey distanceFromHead = {distanceFromHeadString()};
    /// ViewKey for the well transmissibility at this perforation
    dataRepository::ViewKey wellTransmissibility = {wellTransmissibilityString()};
    /// ViewKey for the well transmissibility at this perforation
    dataRepository::ViewKey wellSkinFactor = { wellSkinFactorString() };
  }
  /// ViewKey struct for the Perforation class
  viewKeysPerforation;

protected:
  void postInputInitialization() override;

private:
  /// Linear distance from well head
  real64 m_distanceFromHead;

  /// Well transmissibility at this perforation
  real64 m_wellTransmissibility;

  /// Well skin factor at this perforation
  real64 m_wellSkinFactor;
};

} // namespace geos

#endif // GEOS_MESH_PERFORATION_HPP
