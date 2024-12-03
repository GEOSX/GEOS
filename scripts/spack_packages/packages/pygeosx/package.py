# Copyright 2013-2023 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class Pygeosx(BundlePackage):
    """This is a set of libraries necessary for the pygeosx ATS environment.
    NOTE: It is recommended to download and install python3 and package
    dependencies using your Linux distribution's package manager (Spack
    builds of python tend to be much slower). """

    version('fakeversion')

    depends_on("py-numpy@1.21.0:1.23.4+blas+lapack")
    depends_on('py-mpi4py')
    depends_on('py-virtualenv')
    depends_on('python@3.10:+shared+pic+tkinter+optimizations')
    depends_on("py-scipy")
    depends_on("openblas")
    depends_on("py-matplotlib")
    depends_on("py-sphinx")
    depends_on("py-sphinx-argparse")
    depends_on("py-test")
