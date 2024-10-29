# Copyright 2013-2023 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *

class PyPydataSphinxTheme(PythonPackage):
    """A clean, three-column, Bootstrap-based Sphinx theme by and for the PyData community"""

    homepage = "https://github.com/pydata/pydata-sphinx-theme"
    pypi = "pydata-sphinx-theme/pydata_sphinx_theme-0.13.3.tar.gz"

    version("0.13.3", sha256="827f16b065c4fd97e847c11c108bf632b7f2ff53a3bca3272f63f3f3ff782ecc")

    depends_on("python", type=("build", "run"))
    depends_on("py-setuptools", type="build")