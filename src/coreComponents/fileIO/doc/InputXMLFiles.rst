.. _InputFiles:

###############################################################################
Input Files
###############################################################################


XML
=================================

GEOS is configured via one (or more) `Extensible Markup Language <https://en.wikipedia.org/wiki/XML>`_ (XML) files.
These files contain a set of elements and attributes that closely follow the internal datastructure of GEOS.
When running GEOS, these files are specified using the `-i` argument:

.. code-block:: bash

    geosx -i input.xml



XML Components
------------------------------

The following illustrates some of the key features of a GEOS-format xml file:

.. code-block:: xml

    <?xml version="1.0" ?>

    <Problem>
        <BlockA
            someAttribute="1.234">

            <!-- Some comment -->
            <BlockB
              name="firstNamedBlock"
              anotherAttribute="0"/>
            <BlockB
              name="secondNamedBlock"
              anotherAttribute="1"/>
        </BlockA>
    </Problem>


The two basic components of an xml file are blocks, which are specified using angle brackets ("<BlockA>  </BlockA>"), and attributes that are attached to blocks (attributeName="attributeValue").
Block and attributes can use any ASCII character aside from `<`, `&`, `'`, and `"` (if necessary, use `&lt;`, `&amp;`, `&apos;`, or `&quot;`).
Comments are indicated as follows: `<!-- Some comment -->`.

At the beginning of a GEOS input file, you will find an optional xml declaration (`<?xml version="1.0" ?>`) that is used to indicate the format to certain text editors.
You will also find the root `Problem` block, where the GEOS configuration is placed.
Note that, aside from these elements and commented text, the xml format requires that no other objects exist at the first level.

In the example above, there is a single element within the `Problem` block: `BlockA`.
`BlockA` has an attribute `someAttribute`, which has a value of 1.234, and has three children: a commented string "Some comment" and two instances of `BlockB`.
The `name` attribute is required for blocks that allow multiple instances, and should include a unique string to avoid potential errors.
Where applicable these blocks will be executed in the order in which they are specified in input file.


Input Validation
=================================

The optional `xmlns:xsi` and `xsi:noNamespaceSchemaLocation` attributes in the Problem block can be used to indicate the type of document and the location of the xml schema to the text editor:

.. code-block:: xml

    <Problem
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xsi:noNamespaceSchemaLocation="/path/to/schema.xsd" />

The schema contains a list of xml blocks and attributes that are supported by GEOS, indicates whether a given object is optional or required, and defines the format of the object (string, floating point number, etc.).
A copy of the schema is included in the GEOS source code (/path/to/GEOS/src/coreComponents/schema/schema.xsd).
It can also be generated using GEOS: ``geosx -s schema.xsd``

Many text editors can use the schema to help in the construction of an xml file and to indicate whether it is valid.
Using a validation tool is highly recommended for all users.
The following instructions indicate how to turn on validation for a variety of tools:


xmllint
---------------------------------

xmllint is a command-line tool that is typically pre-installed on UNIX-like systems.
To check whether an input file is valid, run the following command:

xmllint --schema /path/to/schema.xsd input_file.xml


Sublime Text
------------------------------

We recommend using the `Exalt <https://github.com/eerohele/exalt>`_ or `SublimeLinter_xmllint <https://github.com/SublimeLinter/SublimeLinter-xmllint>`_ plug-ins to validate xml files within sublime.
If you have not done so already, install the sublime `Package Control <https://packagecontrol.io/installation>`_.
To install the package, press ``ctrl + shift + p``, type and select ``Package Control: Install Package``, and search for ``exalt`` or ``SublimeLinter`` / ``SublimeLinter-xmllint``.
Note that, depending on the circumstances, these tools may indicate only a subset of the validation errors at a given time.
Once resolved, the tools should re-check the document to look for any additional errors.

As an additional step for SublimLinter-xmllint, you will need to add a linter configuration.
To do so, go to Preferences/Package Settings/SublimeLinter/Settings.
In the right-hand side of the new window, add the xmllint configuration:

.. code-block:: python

    {
        "linters": {
            "xmllint":
            {
                "args": "--schema /path/to/schema.xsd",
                "styles": [
                    {
                        "mark_style": "fill",
                        "scope": "region.bluish",
                        "types": ["error"],
                        "icon": "stop",
                    }
                ]
            },
        }
    }



VS Code
------------------------------

We recommend using the `XML <https://marketplace.visualstudio.com/items?itemName=redhat.vscode-xml>`_ for validating xml files.
After installing this extension, you can associate GEOS format xml files by adding the following entry to the user settings file (replacing `systemId` with the correct path to the schema file):


 .. code-block:: python

    {
        "xml.fileAssociations": [

            {
                "pattern": "**.xml",
                "systemId": "/path/to/GEOS/src/coreComponents/schema/schema.xsd"
            }
        ]
    }


Eclipse
------------------------------

The Eclipse Web Develop Tools includes features for validating xml files.
To install them, go to Help -> Eclipse Marketplace, search for the Eclipse Web Developer Tools, install the package, and restart Eclipse.
Finally, configure the xml validation preferences under Window -> Preferences -> XML -> XML Files -> Validation.
Eclipse will automatically fetch the schema, and validate an active xml file.
The editor will highlight any lines with errors, and underline the specific errors.


GEOS XML Tools
------------------------------

The geosx_xml_tools package, which is used to enable advanced features such as parameters, symbolic math, etc., contains tools for validating xml files.
To do so, call the command-line script with the -s argument, i.e.: `preprocess_xml input_file.xml -s /path/to/schema.xsd`.
After compiling the final xml file, pygeosx will fetch the designated schema, validate, and print any errors to the screen.

Note: Attributes that are using advanced xml features will likely contain characters that are not allowed by their corresponding type pattern.
As such, file editors that are configured to use other validation methods will likely identify errors in the raw input file.


XML Schema
=================================

An XML schema definition (XSD) file lays out the expected structure of an input XML file.
During the build process, GEOS automatically constructs a comprehensive schema from the code's data structure, and updates the version in the source (GEOS/src/coreComponents/schema/schema.xsd).


Schema Components
------------------------------

The first entry in the schema are a set of headers the file type and version.
Following this, the set of available simple types for attributes are laid out.
Each of these includes a variable type name, which mirrors those used in the main code, and a regular expression, which is designed to match valid inputs.
These patterns are defined and documented in ``rtTypes`` (in ``DataTypes.hpp``.
The final part of the schema is the file layout, beginning with the root ``Problem``.
Each complex type defines an element, its children, and its attributes.
Each attribute defines the input name, type, default value, and/or usage.
Comments preceding each attribute are used to relay additional information to the users.


Automatic Schema Generation
------------------------------

A schema may be generated by calling the main code with the -s argument , e.g.: ``geosx -s schema.xsd`` (Note: this is done automatically during the bulid process).
To do this, GEOS does the following:

  1) Initialize the GEOS data structure.
  2) Initialize objects that are registered to catalogs via ``ManagedGroup::ExpandObjectCatalogs()``.
  3) Recursively write element and attribute definitions to the schema using information stored in GEOS groups and wrappers.
  4) Define any expected deviations from the schema via ``ManagedGroup::SetSchemaDeviations()``.


.. _AdvancedXMLFeatures:

Advanced XML Features
=================================

The `geosx_xml_tools` python package adds a set of advanced features to the GEOS xml format: units, parameters, and symbolic expressions.
See`Python Tools Setup <https://geosx-geosx.readthedocs-hosted.com/projects/geosx-geospythonpackages/en/latest/>`_ for details on setup instructions, and `XML Parser Documentation <https://geosx-geosx.readthedocs-hosted.com/projects/geosx-geospythonpackages/en/latest/geosx_xml_tools.html>`_ for package API details.


Usage
---------------------------------

An input file that uses advanced xml features requires preprocessing before it can be used with GEOS.
The preprocessor writes a compiled xml file to the disk, which can be read directly by GEOS and serves as a permanent record for the simulation.
There are three ways to apply the preprocessor:

1) Automatic Preprocessing:  Substituting `geosx` for `geosx_preprocessed` when calling the code will automatically apply the preprocessor to the input xml file, and then pass the remaining arguments to GEOS.  With this method, the compiled xml files will have the suffix '.preprocessed'.  Before running the code, the compiled xml file will also be validated against the xml schema.

.. code-block:: bash

    # Serial example
    geosx_preprocessed -i input.xml

    # Parallel example
    srun -n 2 geosx_preprocessed -i input.xml -x 2


2) Manual Preprocessing:  For this approach, xml files are preprocessed manually by the user with the `preprocess_xml` script.  These files can then be submitted to GEOS separately:

.. code-block:: bash

    # The -c argument is used to manually specify the compiled name
    preprocess_xml -i input.xml -c input.xml.processed
    geosx -i input.xml.processed

    # Otherwise, a random name will be chosen by the tool
    compiled_input=$(preprocess_xml input.xml)
    geosx -i $compiled_input


3) Python / pygeosx: The preprocessor can also be applied directly in python or in pygeosx simulations.  An example of this is method is provided here: `GEOS/examples/pygeosxExamples/hydraulicFractureWithMonitor/`.


Each of these options support specifying multiple input files via the command line (e.g. `geosx_preprocessed -i input_a.xml -i input_b.xml`).
They also support any number of command-line parameter overrides (e.g. `geosx_preprocessed -i input_a.xml -p parameter_a alpha -p parameter_b beta`).


Included Files
------------------------------

Both the XML preprocessor and GEOS executable itself provide the capability to build complex
multi-file input decks by including XML files into other XML files.

The files to be included are listed via the `<Included>` block. There maybe any number of such blocks.
Each block contains a list of `<File name="..."/>` tags, each indicating a file to include.
The `name` attribute must contain either an absolute or a relative path to the included file.
If the path is relative, it is treated as relative to the location of the referring file.
Included files may also contain includes of their own, i.e. it is possible to have `a.xml` include `b.xml`
which in turn includes `c.xml`.

.. note::
   When creating multi-file input decks, it is considered best practice to use relative file paths.
   This applies both to XML includes, and to other types of file references (for example, table file names).
   Relative paths keep input decks both relocatable within the file system and sharable between users.

XML preprocessor's merging capabilities are more advanced than GEOS built-in ones.
Both are outlined below.

XML preprocessor
^^^^^^^^^^^^^^^^

The merging approach is applied recursively, allowing children to include their own files.
Any potential conflicts are handled via the following scheme:

- Merge two objects if:
    - At the root level an object with the matching tag exists.
    - If the "name" attribute is present and an object with the matching tag and name exists.
    - Any preexisting attributes on the object are overwritten by the donor.
- Otherwise append the XML structure with the target.

GEOS
^^^^^

GEOS's built-in processing simply inserts the included files' content (excluding the root node)
into the XML element tree, at the level of `<Included>` tag. Partial merging is handled implicitly
by GEOS's data structure, which treats repeated top-level XML blocks as if they are one single block.
This is usually sufficient for merging together top-level input sections from multiple files,
such as multiple `<FieldSpecifications>` or `<Events>` sections, but more complex cases may require
the use of preprocessor.

.. note::
   While GEOS's XML processing is capable of handling any number of `<Included>` block at any level,
   the XML schema currently produced by GEOS only allows a single such block, and only directly
   within the `<Problem>` tag. Inputs that use multiple blocks or nest them deeper may run but will
   fail to validate against the schema. This is a known discrepancy that may be fixed in the future.

Parameters
------------------------------

Parameters are a convenient way to build a configurable and human-readable input XML.
They are defined via a block in the XML structure.
To avoid conflicts with other advanced features, parameter names can include upper/lower case letters and underscores.
Parameters may have any value, including:

- Numbers (with or without units)
- A path to a file
- A symbolic expression
- Other parameters
- Etc.

They can be used as part of any input xml attribute as follows:

- $x_par$  (preferred)
- $x_par
- $:x_par
- $:x_par$

Attributes can be used across Included files, but cannot be used to set the names of included files themselves.
The following example uses parameters to set the root path for a table function, which is then scaled by another parameter:

.. code-block:: xml

  <Parameters>
    <Parameter
      name="flow_scale"
      value="0.5"/>
    <Parameter
      name="table_root"
      value="/path/to/table/root"/>
  </Parameters>
  
  <FieldSpecifications>
    <SourceFlux
      name="sourceTerm"
      objectPath="ElementRegions/Region1/block1"
      scale="$flow_scale$"
      functionName="flow_rate"
      setNames="{ source }"/>
  </FieldSpecifications>

  <Functions>
    <TableFunction
      name="flow_rate"
      inputVarNames="{time}"
      coordinateFiles="{$table_root$/time_flow.geos}"
      voxelFile="$table_root$/flow.geos"
      interpolation="linear"/>
  </Functions>

Any number of parameter overrides can be issued from the command line using the `-p name value` argument in the preprocessor script.
Note that if the override value contains any spaces, it may need to be surrounded by quotation marks (`-p name "paramter with spaces"`).


Units
------------------------------

The units for any input values to GEOS can be in any self-consistent system.
In many cases, it is useful to override this behavior by explicitly specifying the units of the input.
These are specified by appending a valid number with a unit definition in square braces.
During pre-processing, these units are converted into base-SI units (we plan to support other unit systems in the future).

The unit manager supports most common units and SI prefixes, using both long- and abbreviated names (e.g.: c, centi, k, kilo, etc.).
Units may include predefined composite units (dyne, N, etc.) or may be built up from sub-units using a python syntax (e.g.: [N], [kg*m/s**2]).
Any (or no) amount of whitespace is allowed between the number and the unit bracket.
The following shows a set of parameters with units specified:

.. code-block:: xml

  <Parameters>
    <Parameter name="paramter_a" value="2[m]"/>
    <Parameter name="paramter_b" value="1.2 [cm]"/>
    <Parameter name="paramter_c" value="1.23e4 [bbl/day]"/>
    <Parameter name="paramter_d" value="1.23E-4 [km**2]"/>
  </Parameters>


Please note that the preprocessor currently does not check whether any user-specified units are appropriate for a given input or symbolic expression.


Symbolic Expressions
------------------------------

Input XML files can also include symbolic mathematical expressions.
These are placed within pairs of backticks (\`), and use a limited python syntax.
Please note that parameters and units are evaluated before symbolic expressions.
While symbolic expressions are allowed within parameters, errors may occur if they are used in a way that results in nested symbolic expressions.
Also, note that residual alpha characters (e.g. `sin(`) are removed before evaluation for security.
The following shows an example of symbolic expressions:

.. code-block:: xml

  <Parameters>
    <Parameter name="a" value="2[m]"/>
    <Parameter name="b" value="1.2 [cm]"/>
    <Parameter name="c" value="3"/>
    <Parameter name="d" value="1.23e-4"/>
  </Parameters>
  <Geometry>
    <Box
      name="perf"
      xMin="{`$a$ - 0.2*$b$`, -1e6, -1e6}"
      xMax="{`$c$**2 / $d$`, 1e6, 1e6}" />
  </Geometry>


Validation
------------------------------

Unmatched special characters ($, [, \`, etc.) in the final xml file indicate that parameters, units, or symbolic math were not specified correctly.  
If the prepreprocessor detects these, it will throw an error and exit.
Additional validation of the compiled files can be completed with `preprocess_xml` by supplying the -s argument and the path to the GEOS schema.

