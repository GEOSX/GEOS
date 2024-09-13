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

#ifndef GEOS_PYTHON_PYWRAPPER_HPP_
#define GEOS_PYTHON_PYWRAPPER_HPP_

// Source includes
#include "dataRepository/WrapperBase.hpp"

namespace geos
{
namespace python
{

/**
 *
 */
PyObject * createNewPyWrapper( dataRepository::WrapperBase & wrapper );

/**
 *
 */
PyTypeObject * getPyWrapperType();

} // namespace python
} // namespace geos

#endif
