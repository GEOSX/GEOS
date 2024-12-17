###############################################################################
Code style
###############################################################################

Introduction
============
GEOS is written in standard c++17. In general, target platforms are:

- Linux
- Mac OS X

Currently, our CI/CD system tests on these platforms:

- Ubuntu 18.04, with gcc 8.0 and clang 8.0.0 + cuda10.1.243
- Centos 7.6.1810, with gcc 8.3.1 + cuda10.1.243
- Centos 7.7, with clang 9.0.0
- Mac OS X, with xcode 11.2

Naming Conventions
==================

File Names
----------------
- File names should be `PascalCase <https://en.wikipedia.org/wiki/Camel_case>`__.
- C++ header files are always named with a file extension of  \*.hpp.
- C++ header implementation files, which contain templated or inline function definitions, are always named \*Helpers.hpp.
- C++ source files are always named with a file extension of  \*.cpp.
- C++ class declarations and definitions are contained files with identical names, except for the extensions.
- C++ free function headers and source files are declared/defined in files with identical names, except for the extension.

For example, a class named "Foo" may be declared in a file named "Foo.hpp", with inline/templated functions
defined in "FooHelpers.hpp", with the source implementation contained in Foo.cpp.

.. warning::

  There should not be identical filenames that only differ by case. Some filesystems are not case-sensitive,
  and worse, some filesystems such as MacOSX are case-preserving but not case sensitive.

Function Names
--------------
Function and member function names should be `camelCase <https://en.wikipedia.org/wiki/Camel_case>`__.

Variable Names
--------------
Variables should be `camelCase <https://en.wikipedia.org/wiki/Camel_case>`__.

Member Names
--------------
Member data should be `camelCase <https://en.wikipedia.org/wiki/Camel_case>`__ prefix with "m\_" (i.e. double m_dataVariable;)

Class/Struct Names
------------------
Please use `PascalCase <https://en.wikipedia.org/wiki/Camel_case>`__ for typenames (i.e. classes)

   .. code-block:: c

      class MyClass;


   .. code-block:: c

      class MyClass
      {
        double m_doubleDataMember;
        int m_integerDataMember;
      }

Alias/Typedef Names
-------------------
Alias and typedefs should be the case of the underlying type that they alias. If no clear format is apparent,
as is the case with `double`, then use `camelCase <https://en.wikipedia.org/wiki/Camel_case>`__

Namespace Names
----------------
Namespaces names are all lower `camel case <https://en.wikipedia.org/wiki/Camel_case>`__.

Example
-------
One example of would be a for a class named "Foo", the declaration would be in a header file named "Foo.hpp"

.. code-block:: c

  /*
   * Foo.hpp
   */

  namespace bar
  {

  class Foo
  {
  public:
    Foo();
  private:
    double m_myDouble;
  }
  }

and a source file named "Foo.cpp"

.. code-block:: c

  /*
   * Foo.cpp
   */
  namespace bar
  {
    Foo::Foo():
      m_myDouble(0.0)
    {
      // some constructor stuff
    }
  }

Const Keyword
================

#. All functions and accessors should be declared as "const" functions unless modification to the class is required.
#. In the case of accessors, both a "const" and "non-const" version should be provided.
#. The const keyword should be placed in the location read by the compiler, which is right to left.

The following examples are provided:

   .. code-block:: c

      int a=0; // regular int
      int const b = 0; // const int
      int * const c = &a; // const pointer to non const int
      int const * const d = &b; // const pointer to const int
      int & e = a; // reference to int
      int const & f = b; // reference to const int


Code Format
===========

GEOS applies a variant of the
`BSD/Allman Style <https://en.wikipedia.org/wiki/Indentation_style#Allman_style>`__.
Key points to the GEOS style are:

#. Opening braces (i.e. "{") go on the next line of any control statement, and are not indented from the control statement.
#. NO TABS. Only spaces. In case it isn't clear ... NO TABS!
#. 2-space indentation

   .. code-block:: c

      for( int i=0 ; i<10 ; ++i )
      {
        std::cout << "blah" << std::endl;
      }

#. Try to stay under 100 character line lengths. To achieve this apply these rules in order
#. Align function declaration/definitions/calls on argument list
#. Break up return type and function definition on new line
#. Break up scope resolution operators

   .. code-block:: c

    void
    SolidMechanics_LagrangianFEM::
    TimeStepExplicit( real64 const& time_n,
                      real64 const& dt,
                      const int cycleNumber,
                      DomainPartition * const domain )
    {
      code here
    }

As part of the continuous integration testing, this GEOS code style is enforced via the uncrustify tool.
While quite extensive, uncrustify does not enforce every example of the preferred code style.
In cases where uncrusitfy is unable to enforce code style, it will ignore formatting rules.
In these cases it is acceptable to proceed with pull requests, as there is no logical recourse.

Header Guards
=============
Header guard names should consist of the name `GEOS`, followed by the component name (e.g. dataRepository),
and finally the name of the header file.
All characters in the macro should be capitalized.
