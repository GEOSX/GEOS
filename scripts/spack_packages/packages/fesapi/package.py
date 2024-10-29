# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import sys

from spack.package import *


class Fesapi(CMakePackage):

    url = "https://github.com/F2I-Consulting/fesapi/archive/refs/tags/v2.4.0.0.tar.gz"
    git = "https://github.com/F2I-Consulting/fesapi.git"

    version("2.4.0.0", sha256="a711e8a1218c876a2799f4d05a9820da71eb5503b5d51b834fae98d9fe635381")

    depends_on("hdf5")
    depends_on("boost@1.67.0")
    depends_on("minizip")

    def cmake_args(self):
        spec = self.spec

        cppflags = " ".join(spec.compiler_flags["cppflags"])
        cxxflags = cppflags + " ".join(spec.compiler_flags["cxxflags"])
        cmake_args = [
            self.define('CMAKE_C_COMPILER', spec['mpi'].mpicc),
            self.define('CMAKE_CXX_COMPILER', spec['mpi'].mpicxx),
            self.define('CMAKE_CXX_FLAGS', cxxflags),
            self.define('HDF5_ROOT', spec['hdf5'].prefix),
            # fesAPI/spack can detect wrong version otherwise
            self.define('HDF5_VERSION', spec['hdf5'].version),
            self.define('MINIZIP_INCLUDE_DIR', spec['minizip'].prefix.include + "/minizip"),
            self.define('MINIZIP_LIBRARY_RELEASE', spec['minizip'].prefix.lib),
            self.define('Boost_INCLUDE_DIR', spec['boost'].prefix.include),
            "-DWITH_EXAMPLE:BOOL=OFF",
            "-DWITH_DOTNET_WRAPPING:BOOL=OFF",
            "-DWITH_JAVA_WRAPPING:BOOL=OFF",
            "-DWITH_PYTHON_WRAPPING:BOOL=OFF",
            "-DWITH_RESQML2_2:BOOL=OFF",
            "-DWITH_TEST:BOOL=OFF",
        ]

        return cmake_args
