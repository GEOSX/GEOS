########################
GEOS Documentation
########################

GEOS is a code framework focused on enabling streamlined development of
physics simulations on high performance computing platforms.  Our documentation
is organized into several separate guides, given that different users will have
different needs.

We recommend all users begin with the Quick Start guide, which
will walk you through downloading and compiling the code.  Application focused
users may then want to explore our Tutorials, which provide an introduction to the
basic capabilities of the code.  More detailed descriptions of these capabilities can
then be found in the User Guide.

For those interested in developing new capabilities in GEOS, we provide a Developer Guide.
The code itself is also documented inline using doxygen.  The Build Guide
contains more detailed information about third-party dependencies, the build
system, and the continuous integration system.  Finally, GEOS has a self-documenting
data structure.  The Datastructure Index is an automatically generated list of
all available input parameters and data structures in the code.  This is a
comprehensive resource, but probably not the place to start.

High quality documentation is a critical component of a successful code.  If
you have suggestions for improving the guides below, please post an issue on our
`issue tracker <https://github.com/GEOS-DEV/GEOS/issues>`_.



.. grid:: 2
    :gutter: 4

    .. grid-item-card::

        Quick Start Guide
        ^^^^^^^^^^^^^^^^^^^

        New to GEOS?  We will walk you through downloading the source, compiling the code, and testing the installation.

        +++

        .. button-ref:: QuickStart
            :expand:
            :color: info
            :click-parent:

            To the Quick Start

    .. grid-item-card::

        Tutorials
        ^^^^^^^^^^

        Working tutorials that show how to run some common problems. After going through these examples, you should have a good understanding of how to set up and solve your own models.

        +++

        .. button-ref:: Tutorials
            :expand:
            :color: info
            :click-parent:

            To the Tutorials

    .. grid-item-card::

        Basic Examples
        ^^^^^^^^^^^^^^^

        Example problems that are organized around physical processes (fluid flow, mechanics, etc.).

        +++

        .. button-ref:: BasicExamples
            :expand:
            :color: info
            :click-parent:

            To the Basic Examples

    .. grid-item-card::

        Advanced Examples
        ^^^^^^^^^^^^^^^^^^^

        Example problems that demonstrate additional physical models, constitutive models, advanced features, etc.

        +++

        .. button-ref:: AdvancedExamples
            :expand:
            :color: info
            :click-parent:

            To the Advanced Examples

    .. grid-item-card::

        User Guide
        ^^^^^^^^^^^^^^^^^^^

        Detailed instructions on how to construct input files, configure problems, manage outputs, etc.

        +++

        .. button-ref:: UserGuide
            :expand:
            :color: info
            :click-parent:

            To the User Guide

    .. grid-item-card::

        Python Tools
        ^^^^^^^^^^^^^^^^^^^

        Documentation for the python packages distributed alongside GEOS used to manage xml files, condition numerical meshes, read outputs, etc.

        +++

        .. button-link:: https://geosx-geosx.readthedocs-hosted.com/projects/geosx-geospythonpackages/en/latest/
            :expand:
            :color: info
            :click-parent:

            To the Python Tools Documentation

    .. grid-item-card::

        Feature Requests, Reporting Bugs, and Support
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        To make feature requests, report bugs, or get support (after reviewing the user guide) please submit an issue on Github.

        +++

        .. button-link:: https://github.com/GEOS-DEV/GEOS/issues/new/choose
            :expand:
            :color: info
            :click-parent:

            To the "New issue" page on the GEOS Github repository


********************
Table of Contents
********************

.. toctree::
   :maxdepth: 2

   docs/sphinx/QuickStart

   docs/sphinx/tutorials/Index

   docs/sphinx/basicExamples/Index

   docs/sphinx/advancedExamples/Index
   
   docs/sphinx/userGuide/Index

   docs/sphinx/developerGuide/Index

   docs/sphinx/Doxygen

   docs/sphinx/buildGuide/Index

   docs/sphinx/CompleteXMLSchema

   docs/sphinx/Contributors

   docs/sphinx/Publications

   docs/sphinx/Acknowledgements

*********************
Indices and tables
*********************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
