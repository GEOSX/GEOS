# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import shutil

# Add python modules to be documented
python_root = './coreComponents/python/modules'
python_modules = ('geos-mesh-tools',
                  'geos-xml-tools',
                  'hdf5-wrapper',
                  'pygeos-tools',
                  'geos-timehistory')
for m in python_modules:
    sys.path.insert(0, os.path.abspath(os.path.join(python_root, m)))

# Call doxygen in ReadtheDocs
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
if read_the_docs_build:

    # Make sure directory exists
    cwd = os.getcwd()

    build_path = os.path.join(cwd, "../_readthedocs")
    if not os.path.isdir(build_path):
        os.mkdir(build_path)

    html_path = os.path.join(build_path, "html")
    if not os.path.isdir(html_path):
        os.mkdir(html_path)

    docs_path = os.path.join(cwd, "docs", "doxygen")
    common_path = os.path.join(cwd, "coreComponents", "common")

    doxyfile_src = os.path.join(docs_path, "Doxyfile.in")
    doxyfile_dst = os.path.join(build_path, "Doxyfile")
    config_src = os.path.join(docs_path, "GeosxConfig.hpp")
    config_dst = os.path.join(common_path, "GeosxConfig.hpp")

    input_dirs = [
        "coreComponents/common",
        "coreComponents/dataRepository",
        "coreComponents/fileIO",
        "coreComponents/linearAlgebra",
        "coreComponents/mesh",
        "coreComponents/finiteElement/elementFormulations",
        "coreComponents/finiteElement/kernelInterface",
        "coreComponents/mesh/MeshFields.hpp",
        "coreComponents/physicsSolvers",
        "coreComponents/finiteVolume",
        "coreComponents/functions",
        "coreComponents/fieldSpecification",
        "coreComponents/discretizationMethods",
        "coreComponents/events",
        "coreComponents/mainInterface"
        ]
        

    # Write correct ReadtheDocs path and input directories
    shutil.copy(doxyfile_src, doxyfile_dst)
    with open(doxyfile_dst, "a") as f:
        f.write("\nINPUT = %s" % " ".join(input_dirs))
        f.write("\nOUTPUT_DIRECTORY = %s/doxygen_output" % html_path)
        f.write("\nHAVE_DOT = YES")

    # Make a symlink to GeosxConfig.hpp in common
    if not os.path.exists(config_dst):
        os.symlink(config_src, config_dst)

    print("********** Running Doxygen in ReadtheDocs **********")
    # Call doxygen
    from subprocess import run
    run(['doxygen', doxyfile_dst])
    print("********** Finished Running Doxygen in ReadtheDocs **********")


# -- Project information -----------------------------------------------------

project = u'GEOS'
copyright = u'2016-2024 Lawrence Livermore National Security LLC, 2018-2024 Total Energies, The Board of Trustees of the Leland Stanford Junior University, 2023-2024 Chevron, 2019- GEOS/GEOSX Contributors'
author = u'GEOS/GEOSX Contributors'

# The short X.Y version
version = u''
# The full version, including alpha/beta/rc tags
release = u''


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_design',
    'sphinx.ext.todo',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.imgmath',
    'sphinxarg.ext',
    'matplotlib.sphinxext.plot_directive',
    'sphinx.ext.napoleon',
    'sphinxcontrib.plantuml',
    'sphinxcontrib.programoutput'
]

plantuml = "/usr/bin/java -Djava.awt.headless=true -jar /tmp/plantuml.jar"
plantuml_output_format = "svg_img"

plot_html_show_source_link = True
plot_html_show_formats = False

autodoc_mock_imports = ["pygeosx", "pylvarray", "meshio", "lxml", "mpi4py", "h5py"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = [u'_build', 'Thumbs.db', '.DS_Store', 'cmake/*', '**/blt/**']

todo_include_todos = True

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Theme options ----------------------------------------------
extensions += [
    'sphinx_rtd_theme',
]

html_theme = "sphinx_rtd_theme"
# html_theme = "pydata_sphinx_theme"

html_theme_options = {
    'navigation_depth': -1,
    'collapse_navigation': False
}


# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
#html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

html_static_path = ['./docs/sphinx/_static']

html_css_files = [
    'theme_overrides.css',
]


# -- Options for HTMLHelp output ---------------------------------------------
# Output file base name for HTML help builder.
htmlhelp_basename = 'GEOSXdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'GEOS.tex', u'GEOS Documentation',
     u'GEOS/GEOSX Developers', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'geos', u'GEOS Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'GEOS', u'GEOS Documentation',
     author, 'GEOS', 'GEOS simulation framework.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------

# Added to allow figure numbering
numfig = True

# Additional stuff for the LaTeX preamble.
latex_elements['preamble'] = '\\usepackage{amsmath}\n\\usepackage{amssymb}\n\\usepackage[retainorgcmds]{IEEEtrantools}\n'


#####################################################
# add LaTeX macros

f = open('docs/sphinx/latex_macros.sty')
imgmath_latex_preamble = ""
imgmath_image_format = 'svg'
imgmath_font_size = 14

for macro in f:
    # used when building latex and pdf versions
    latex_elements['preamble'] += macro + '\n'
    # used when building html version
    imgmath_latex_preamble += macro + '\n'

#####################################################
