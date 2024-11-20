# Copyright 2013-2023 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *

class PySphinxcontribPlantuml(PythonPackage):
    """PlantUML for Sphinx"""

    homepage = "https://github.com/sphinx-contrib/plantuml"
    pypi = "sphinxcontrib-plantuml/sphinxcontrib-plantuml-0.26.tar.gz"

    version("0.26", sha256="adb3397d5cb0613632cd3dad7894381422bac24464c393cb050404dd6712b1a7")

    depends_on("python", type=("build", "run"))
    depends_on("py-setuptools", type="build")