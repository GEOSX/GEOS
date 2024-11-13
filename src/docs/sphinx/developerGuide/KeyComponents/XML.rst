.. _XML_and_classes:

XML Input
##################################################


In this document, you will learn how GEOS classes interact with external information parsed from XML files, and how to add a new XML block that can be interpreted by GEOS.
Flow solvers and relative permeability are used as examples.


GEOS data structure overview
=============================

.. _GroupPar:

Group : the base class of GEOS
-------------------------------

All GEOS classes derive from a base class called ``dataRepository::Group``.
The ``Group`` class provides a way to organize all GEOS objects in a filesystem-like structure.
One could think of ``Group`` s as *file folders* that can bear data (stored in ``Wrapper`` s), have a parent folder (another ``Group``),  and have possibly multiple subfolders (referred to as the subgroups).
Below, we briefly review the data members of the ``Group`` class that are essential to understand the correspondence between the GEOS data structure and the XML input.
For more details, we refer the reader to the extensive documentation of the :ref:`dataRepository`, including the :ref:`Group` class documentation.


In the code listing below, we see that each ``Group`` object is at minimum equipped with the following member properties:

- A pointer to the parent ``Group`` called ``m_parent`` (member classes are prefixed by ``m_``),
- The ``Group`` 's own data, stored for flexibility in an array of generic data ``Wrapper`` s called ``m_wrappers``,
- A map of one or many children (also of type ``Group``) called ``m_subGroups``.
- The ``m_size`` and ``m_capacity`` members, that are used to set the size and capacity of any objects contained.
- The name of the ``Group``, stored as a ``string`` in ``m_name``. This name can be seen as the object unique ID.

.. literalinclude:: ../../../../coreComponents/dataRepository/Group.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_02
   :end-before: RestartFlags

*[Source: src/coreComponents/dataRepository/Group.hpp]*

.. _ObjectCatalogPar:

A few words about the ObjectCatalog
-----------------------------------

**What is an ObjectCatalog and why do we need it?**

Some classes need external information (physical and/or algorithmic parameters for instance) provided by the user to be instantiated.
This is the case when the ``m_input_flags`` data member of one of the ``Group`` 's ``Wrapper`` s has an entry set to ``REQUIRED`` (we will illustrate this below).
In this situation, the required information must be supplied in the XML input file, and if it is absent, an error is raised by GEOS.

To connect the external (XML) and internal (C++) data structures, GEOS uses an **ObjectCatalog** that maps keys (of type ``string``) to the corresponding classes (one unique key per mapped class).
These string keys, referred to as ``catalogName`` s, are essential to transfer the information from the XML file to the factory functions in charge of object instantiation (see below).

**What is a CatalogName?**

The ``catalogName`` of an object is a *key* (of type ``string``) associated with this object's class.
On the one hand, in the XML file, the key is employed by the user as an XML tag to specify the type of object (e.g., the type of solver, constitutive model, etc) to create and use during the simulation.
On the other hand, internally, the key provides a way to access the appropriate factory function to instantiate an object of the desired class.

Most of the time, the ``catalogName`` and the C++ class name are identical.
This helps make the code easier to debug and allows the XML/C++ correspondence to be evident.
But strictly speaking, the ``catalogName`` can be anything, as long as it refers uniquely to a specific class.
The ``catalogName`` must not be confused with the object's *name* (``m_name`` is a data member of the class that stores the object's unique ID, not its class key).
You can have several objects of the same class and hence the same ``catalogName``, but with different names (i.e. unique ID): several fluid models, several solvers, etc.

**How can I add my new externally-accessible class to the ObjectCatalog?**

Let us consider a flow solver class derived from ``FlowSolverBase``, that itself is derived from ``PhysicsSolverBase``.
To instantiate and use this solver, the developer needs to make the derived flow solver class reachable from the XML file, via an XML tag.
Internally, this requires adding the derived class information to ``ObjectCatalog``, which is achieved with two main ingredients: 1) a ``CatalogName()`` method in the class that lets GEOS know *what* to search for in the internal ``ObjectCatalog`` to instantiate an object of this class, 2) a macro that specifies *where* to search in the ``ObjectCatalog``.

1. To let GEOS know what to search for in the catalog to instantiate an object of the derived class, the developer must equip the class  with a ``CatalogName()`` method that returns a ``string``.
   In this document, we have referred to this returned ``string`` as the object's ``catalogName``, but in fact, the method ``CatalogName()`` is what matters since the ``ObjectCatalog`` contains all the ``CatalogName()`` return values.
   Below, we illustrate this with the ``CompositionalMultiphaseFlow`` solver.
   The first code listing defines the class name, which in this case is the same as the ``catalogName`` shown in the second listing.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: public:

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: virtual

*[Source: src/coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp]*


2. To let GEOS know where to search in the ``ObjectCatalog``, a macro needs to be added at the end of the .cpp file implementing the class.
   This macro (illustrated below) must contain the type of the base class (in this case, ``PhysicsSolverBase``), and the name of the derived class (continuing with the example used above, this is ``CompositionalMultiphaseFlow``).
   As a result of this construct, the ``ObjectCatalog`` is not a flat list of ``string`` s mapping the C++ classes.
   Instead, the ``ObjectCatalog`` forms a tree that reproduces locally the structure of the class diagram, from the base class to the derived classes.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: //END_SPHINX_INCLUDE_01

*[Source: src/coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.cpp]*


.. highlights::
  Summary: All GEOS objects form a filesystem-like structure. If an object needs to be accessible externally, it must be registered in the ``ObjectCatalog``. This is done by adding ``CatalogName()`` method that returns a ``string`` key to the object's class, and by adding the appropriate macro. The catalog has the same tree structure as the class diagram.


Registration: parsing XML input files to instantiate GEOS objects
------------------------------------------------------------------


In this section, we describe with more details the connection between **internal GEOS objects** and **external XML tags** parsed from parameter files.
We call this process *Registration*.
The registration process works in three steps:

  #. The XML document is parsed.
     Each time a new XML tag is found, the current local scope of the ``ObjectCatalog`` is inspected.
     The goal is to find a ``catalogName`` ``string`` that matches the XML tag.
  #. If it is the case (the current local scope of the ``ObjectCatalog`` contains a ``catalogName`` identical to the XML tag), then the code creates a new instance of the class that the ``catalogName`` refers to.
     This new object is inserted in the ``Group`` tree structure at the appropriate location, as a subgroup.
  #. By parsing the XML attributes of the tag, the new object properties are populated.
     Some checks are performed to ensure that the data supplied is conform, and that all the required information is present.


Let's look at this process in more details.

Creating a new object and giving it a Catalog name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider again that we are registering a flow solver deriving from ``FlowSolverBase``, and assume that this solver is called ``CppNameOfMySolver``.
This choice of name is not recommended (we want names that reflect what the solver does!), but for this particular example, we just need to know that this name is the class name inside the C++ code.

To specify parameters of this new solver from an XML file, we need to be sure that the XML tag and the ``catalogName`` of the class are identical.
Therefore, we equip the ``CppNameOfMySolver`` class with a ``CatalogName()`` method that returns the solver ``catalogName`` (=XML name).
Here, this method returns the ``string``  "XmlNameOfMySolver".

We have deliberately distinguished the class name from the catalog/XML name for the sake of clarity in this example.
It is nevertheless a best practice to use the same name for the class and for the ``catalogName``.
This is the case below for the existing ``CompositionalMultiphaseFVM`` class.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: //END_SPHINX_INCLUDE_01

*[Source: src/coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp]*


Parsing XML and searching the ObjectCatalog in scope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that we have implemented a ``CatalogName()`` method returning a specific key (of type ``string``), we can have a block in our XML input file with a tag that corresponds to the ``catalogName`` "XmlNameOfMySolver".
This is how the XML block would look like.

.. code-block:: xml

    <Problem>
      <Solvers
        gravityVector="{ 0.0, 0.0, -9.81 }">
        <XmlNameOfMySolver name="nameOfThisSolverInstance"
                                 verboseLevel="1"
                                 gravityFlag="1"
                                 temperature="297.15" />
          <LinearSolverParameters newtonTol="1.0e-6"
                                  maxIterNewton="15"
                                  useDirectSolver="1"/>
        </XmlNameOfMySolver>
      </Solvers>
    </Problem>

Here, we see that the XML structure defines a parent node "Problem", that has (among many others) a child node "Solvers".
In the "Solvers" block, we have placed the new solver block as a child node of the "Solvers" block with the XML tag corresponding to the ``catalogName`` of the new class.
We will see in details next how the GEOS internal structure constructed from this block mirrors the XML file structure.


Instantiating the new solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Above, we have specified an XML block with the tag "XmlNameOfMySolver".
Now, when reading the XML file and encountering an "XmlNameOfMySolver" solver block, we add a new instance of the class ``CppNameOfMySolver`` in the filesystem structure as explained below.

We saw that in the XML file, the new solver block appeared as child node of the XML block "Solvers".
The internal construction mirrors this XML structure.
Specifically, the new object of class ``CppNameOfMySolver`` is registered as a subgroup (to continue the analogy used so far, as a subfolder) of its parent ``Group``, the class ``PhysicsSolverManager`` (that has a ``catalogName`` "Solvers").
To do this, the method ``CreateChild`` of the ``PhysicsSolverManager`` class is used.


.. code-block:: cpp

    // Variable values in this example:
    // --------------------------------
    // childKey = "XmlNameOfMySolver" (string)
    // childName = "nameOfThisSolverInstance" (string)
    // PhysicsSolverBase::CatalogInterface = the Catalog attached to the base Solver class
    // hasKeyName = bool method to test if the childKey string is present in the Catalog
    // registerGroup = method to create a new instance of the solver and add it to the group tree

.. literalinclude:: ../../../../coreComponents/physicsSolvers/PhysicsSolverManager.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: void

*[Source: src/coreComponents/physicsSolvers/PhysicsSolverManager.cpp]*

In the code listing above, we see that in the ``PhysicsSolverManager`` class, the ``ObjectCatalog`` is searched to find the ``catalogName`` "CompositionalMultiphaseFlow" in the scope of the ``PhysicsSolverBase`` class.
Then, the factory function of the base class ``PhysicsSolverBase`` is called.
The ``catalogName`` (stored in ``childKey``) is passed as an argument of the factory function to ensure that it instantiates an object of the desired derived class.

As explained above, this is working because 1) the XML tag matches the ``catalogName`` of the ``CompositionalMultiphaseFlow`` class and 2) a macro is placed at the end of the .cpp file implementing the ``CompositionalMultiphaseFlow`` class to let the ``ObjectCatalog`` know that ``CompositionalMultiphaseFlow`` is a derived class of ``PhysicsSolverBase``.

Note that several instances of the same type of solver can be created, as long as they each have a different name.


Filling the objects with data (wrappers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After finding and placing the new solver ``Group`` in the filesystem hierarchy, properties are read and stored.
This is done by registering *data wrappers*.
We refer to the documentation of the :ref:`dataRepository` for additional details about the ``Wrapper`` s.
The method used to do that is called ``registerWrapper`` and is placed in the class constructor when the data is required in the XML file.
Note that some properties are registered at the current (derived) class level, and other properties can also be registered at a base class level.

Here, the only data (=wrapper) that is defined at the level of our ``CppNameOfMySolver`` class is temperature, and everything else is registered at the base class level.
We register a property of temperature, corresponding to the member class ``m_temperature`` of ``CppNameOfMySolver``.
The registration also checks if a property is required or optional (here, it is required), and provides a brief description that will be used in the auto-generated code documentation.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseBase.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: Mass

*[Source: src/coreComponents/physicsSolvers/fluidFlow/CompositionalMultiphaseBase.cpp]*

This operation is done recursively if XML tags are nested.


To summarize:
-------------


  - Every class in GEOS derive from a ``Group`` in a filesystem-like structure.
    A ``Group`` must have a parent ``Group``, can have data (in ``Wrapper`` s), and can have one or many children (the subgroups).
    There is an ``ObjectCatalog`` in which the classes derived from ``Group`` are identified by a key called the ``catalogName``.
  - When parsing XML input files, GEOS inspects each object's scope in the ``ObjectCatalog`` to find classes with the same ``catalogName`` as the XML tag.
    Once it finds an XML tag in the ``ObjectCatalog``, it registers it inside the filesystem structure.
  - In the registration process, properties from the XML file are parsed and used to allocate member data ``Wrapper`` s and fully instantiate the ``Group`` class.
  - If XML tags are nested, subgroups are allocated and processed in a nested manner.

The correspondence between XML and class hierarchy is thus respected, and the internal object hierarchy mirrors the XML structure.



Example: adding a new relative permeability model
=================================================

This example is taken from the class ``BrooksCoreyRelativePermeability``, derived from ``RelativePermeabilityBase``.


Implement a ``CatalogName`` function (.hpp):
--------------------------------------------

As explained above we add the class to the ``ObjectCatalog`` in two steps. First we implement the ``CatalogName`` function:

.. literalinclude:: ../../../../coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: virtual

*[source: src/coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.hpp]*

Then in the .cpp file we add the macro to register the catalog entry:

.. literalinclude:: ../../../../coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: }

*[source: src/coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.cpp]*

Now every time a "BrooksCoreyRelativePermeability" ``string`` is encountered inside a ``Relative Permeability`` catalog, we will instantiate a class ``BrooksCoreyRelativePermeability``.

Declare the ``Wrapper`` s  keys (.hpp):
---------------------------------------

When attaching properties (i.e. data ``Wrapper`` s) to a class, a similar registration process must be done.
Every property is accessed through its ``ViewKey`` namespace.
In this namespace, we define ``string`` s that correspond to the tags of XML attributes of the "BrooksCoreyRelativePermeability" block.

.. literalinclude:: ../../../../coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: //END_SPHINX_INCLUDE_01

*[source: src/coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.hpp]*

Declare data members (.hpp):
----------------------------

The data members are defined in the class.
They will ultimately contain the data read from the XML file (other data members not read from the XML file can also exist).

.. literalinclude:: ../../../../coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_02
   :end-before: //END_SPHINX_INCLUDE_02

*[source: src/coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.hpp]*

Implement the data registration process (``registerWrapper``):
--------------------------------------------------------------

The registration process done in the class constructor puts everything together.
It connects the attributes values in the XML file to class member data.
For instance, in the listing below, the first ``registerWrapper`` call means that we want to read in the XML file the attribute value corresponding to the attribute tag ''phaseMinVolumeFraction'' defined in the .hpp file, and that we want to store the read values into the ``m_phaseMinVolumeFraction`` data members.
We see that this input is not required.
If it is absent from the XML file, the default value is used instead.
The short description that completes the registration will be added to the auto-generated documentation.

.. literalinclude:: ../../../../coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: }

*[source: src/coreComponents/constitutive/relativePermeability/BrooksCoreyRelativePermeability.cpp]*

The XML block
-------------

We are ready to use the relative permeability model in GEOS.
The corresponding XML block (child node of the "Constitutive" block) reads:

.. code-block:: XML

  <Constitutive>
    <BrooksCoreyBakerRelativePermeability name="relperm"
                                          phaseNames="{oil, gas, water}"
                                          phaseMinVolumeFraction="{0.05, 0.05, 0.05}"
                                          waterOilRelPermExponent="{2.5, 1.5}"
                                          waterOilRelPermMaxValue="{0.8, 0.9}"
                                          gasOilRelPermExponent="{3, 3}"
                                          gasOilRelPermMaxValue="{0.4, 0.9}"/>
  <Constitutive>

With this construct, we instruct the ``ConstitutiveManager`` class (whose ``catalogName`` is "Constitutive") to instantiate a subgroup of type ``BrooksCoreyRelativePermeability``.
We also fill the data members of the values that we want to use for the simulation.
For a simulation with multiple regions, we could define multiple relative permeability models in the "Constitutive" XML block (yielding multiple relperm subgroups in GEOS), with a unique name attribute for each model.

*For more examples on how to contribute to GEOS, please read* :ref:`AddingNewSolver`

.. _`InputSchemaGenerationPar`:

Input Schema Generation
===========================

A schema file is a useful tool for validating input .xml files and constructing user-interfaces.  Rather than manually maintaining the schema during development, GEOS is designed to automatically generate one by traversing the documentation structure.

To generate the schema, run GEOS with the input, schema, and the (optional) schema_level arguments, i.e.: ``geosx -i input.xml -s schema.xsd``.  There are two ways to limit the scope of the schema:

1. Setting the verbosity flag for an object in the documentation structure.  If the schema-level argument is used, then only objects (and their children) and attributes with ``(verbosity < schema-level)`` will be output.

2. By supplying a limited input xml file.  When GEOS builds its data structure, it will only include objects that are listed within the xml (or those that are explicitly appended when those objects are initialized).  The code will add all available *attributes* for these objects to the schema.

To take advantage of this design it is necessary to use the automatic xml parsing approach that relies upon the documentation node.  If values are read in manually, then the schema can not be used to validate xml those inputs.

Note: the lightweight xml parser that is used in GEOS cannot be used to validate inputs with the schema directly.  As such, it is necessary to use an external tool for validation, such as the geosx_tools python module.
