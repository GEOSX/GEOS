# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import sys

from spack import *


class Vtk(CMakePackage):
    """The Visualization Toolkit (VTK) is an open-source, freely
    available software system for 3D computer graphics, image
    processing and visualization. """

    homepage = "http://www.vtk.org"
    url = "https://www.vtk.org/files/release/9.0/VTK-9.0.0.tar.gz"
    list_url = "http://www.vtk.org/download/"

    maintainers = ['chuckatkins', 'danlipsa']

    version("9.2.6", sha256="06fc8d49c4e56f498c40fcb38a563ed8d4ec31358d0101e8988f0bb4d539dd12")
    version('9.1.0', sha256='8fed42f4f8f1eb8083107b68eaa9ad71da07110161a3116ad807f43e5ca5ce96')
    version('9.0.3', sha256='bc3eb9625b2b8dbfecb6052a2ab091fc91405de4333b0ec68f3323815154ed8a')
    version('9.0.1', sha256='1b39a5e191c282861e7af4101eaa8585969a2de05f5646c9199a161213a622c7')
    version('9.0.0', sha256='15def4e6f84d72f82386617fe595ec124dda3cbd13ea19a0dcd91583197d8715')
    version('8.2.0', sha256='34c3dc775261be5e45a8049155f7228b6bd668106c72a3c435d95730d17d57bb')
    version('8.1.2', sha256='0995fb36857dd76ccfb8bb07350c214d9f9099e80b1e66b4a8909311f24ff0db')
    version('8.1.1', sha256='71a09b4340f0a9c58559fe946dc745ab68a866cf20636a41d97b6046cb736324')
    version('8.1.0', sha256='6e269f07b64fb13774f5925161fb4e1f379f4e6a0131c8408c555f6b58ef3cb7')
    version('8.0.1', sha256='49107352923dea6de05a7b4c3906aaf98ef39c91ad81c383136e768dcf304069')
    version('7.1.0', sha256='5f3ea001204d4f714be972a810a62c0f2277fbb9d8d2f8df39562988ca37497a')
    version('7.0.0', sha256='78a990a15ead79cdc752e86b83cfab7dbf5b7ef51ba409db02570dbdd9ec32c3')
    version('6.3.0', sha256='92a493354c5fa66bea73b5fc014154af5d9f3f6cee8d20a826f4cd5d4b0e8a5e')
    version('6.1.0', sha256='bd7df10a479606d529a8b71f466c44a2bdd11fd534c62ce0aa44fad91883fa34')

    variant('python', default=False, description='Enable Python support')
    variant('mpi', default=True, description='Enable MPI support')

    patch('vtkXMLReader-fpe.patch', when='@9.1.0:')

    extends('python', when='+python')

    # Acceptable python versions depend on vtk version
    # We need vtk at least 8.0.1 for python@3,
    # and at least 9.0 for python@3.8
    depends_on('python@2.7:2.9', when='@:8.0 +python', type=('build', 'run'))
    depends_on('python@2.7:3.7.99', when='@8.0.1:8.9 +python', type=('build', 'run'))
    depends_on('python@2.7:', when='@9.0: +python', type=('build', 'run'))

    # We need mpi4py if buidling python wrappers and using MPI
    depends_on('py-mpi4py', when='+python+mpi', type='run')

    depends_on('mpi', when='+mpi')

    def cmake_args(self):
        spec = self.spec

        # yapf: disable
        # Added GEOSX Arguments
        if '+mpi' in spec:
            mpi_args = [
                self.define('CMAKE_C_COMPILER', spec['mpi'].mpicc),
                self.define('CMAKE_CXX_COMPILER', spec['mpi'].mpicxx),
                self.define('CMAKE_CXX_FLAGS', self.spec.compiler_flags["cxxflags"]),
                '-DVTK_USE_MPI=ON',
                '-DVTK_MODULE_ENABLE_VTK_IOParallelXML=YES',
                '-DVTK_MODULE_ENABLE_VTK_FiltersParallelDIY2=YES'
            ]
        else:
            mpi_args = [
                self.define('CMAKE_C_COMPILER', self.compiler.cc),
                self.define('CMAKE_CXX_COMPILER', self.compiler.cxx),
                self.define('CMAKE_CXX_FLAGS', self.spec.compiler_flags["cxxflags"]),
                '-DVTK_USE_MPI=OFF',
                '-DVTK_MODULE_ENABLE_VTK_IOParallelXML=NO',
                '-DVTK_MODULE_ENABLE_VTK_FiltersParallelDIY2=NO',
            ]

        cmake_args= [
            '-DVTK_GROUP_ENABLE_Imaging=DONT_WANT',
            '-DVTK_GROUP_ENABLE_MPI=DONT_WANT',
            '-DVTK_GROUP_ENABLE_Qt=DONT_WANT',
            '-DVTK_GROUP_ENABLE_Rendering=DONT_WANT',
            '-DVTK_GROUP_ENABLE_StandAlone=DONT_WANT',
            '-DVTK_GROUP_ENABLE_Views=DONT_WANT',
            '-DVTK_GROUP_ENABLE_Web=DONT_WANT',
            '-DVTK_BUILD_ALL_MODULES=OFF',
            '-DVTK_WRAP_PYTHON=OFF',
            '-DVTK_WRAP_JAVA=OFF',
            '-DVTK_MODULE_ENABLE_VTK_vtkm=DONT_WANT',
            '-DVTK_MODULE_ENABLE_VTK_IOXML=YES',
            '-DVTK_MODULE_ENABLE_VTK_IOLegacy=YES',
            '-DVTK_BUILD_TESTING=OFF',
            '-DVTK_LEGACY_REMOVE=ON'
        ]
        # yapf: enable

        cmake_args = mpi_args + cmake_args

        return cmake_args
