.. _Doxygen:

Doxygen Documentation
##################################################

Developer documentation of code is provided in the form of Doxygen-style comment blocks.
`Doxygen <http://www.doxygen.nl/>`__ is a tool for generating html/xml/latex documentation for C++ code from specially marked code comments.
Having concise but high quality documentation of public APIs helps both users and developers of these APIs.
We use Doxygen and Ctest to enforce documentation coverage. If Doxygen produces any warnings, your pull request will fail CI checks!
See :ref:`GitWorkflow` for more on pull requests and CI.

Accessing
====================================

There are two ways to access Doxygen documentation.

Build locally
------------------------------------

Prior to configuring a GEOS build, have Doxygen installed:

  .. code-block:: sh

   sudo apt install doxygen

.. note::

  Eventually, doxygen (version 1.8.13) is provided within the `thirdPartyLibs` repository.

Configure GEOS and go the build directory:

  .. code-block:: sh

   cd GEOS/build-your-platform-release

Build doxygen docs only:

  .. code-block:: sh

   make geosx_doxygen

Or build all docs:

  .. code-block:: sh

   make geosx_docs

Open in browser:

  .. code-block:: sh

   google-chrome html/doxygen_output/html/index.html

On readthedocs
------------------------------------

Go to `GEOS documentation <https://geosx-geosx.readthedocs-hosted.com/>`__, select the version of interest, and follow
the Doxygen link at the left-hand-side.

Guidelines
====================================

What to document
------------------------------------

The following entities declared in project header files within `geosx` namespace require documentation:

- all classes and structs, including public nested ones
- global functions, variables and type aliases
- public and protected member functions, variables and type aliases in classes
- preprocessor macros

Exceptions are made for:

- overrides of virtual functions in derived types
- implementation details nested in namespace `internal`
- template specializations in some cases

How to document
------------------------------------

The following rules and conventions are used. Some are stricter than others.

#. We use `@`-syntax for all Doxygen commands (e.g. `@brief` instead of `\\brief`).

#. Entities such as type aliases and member variables that typically only require a brief description,
   can have a single-line documentation starting with `///`.

   * `@brief` is not required for single-line comments.

#. Entities such as classes and functions that typically require either detailed explanation or parameter documentation,
   are documented with multiline comment blocks.

   * `@brief` is required for comment blocks.

#. Brief and detailed descriptions should be complete sentences (i.e. start with a capital letter and end with a dot).

#. Prefer concise wording in `@brief`, e.g. "Does X." instead of "This is a function that does X."

#. All functions parameters and return values must be explicitly documented via `@param` and `@return`.

   * An exception to this rule seem to be copy/move constructor/assignment, where parameter documentation can be omitted.

#. Add `[in]` and `[out]` tags to function parameters, as appropriate.

#. Function and template parameter descriptions are not full sentences (i.e. not capitalized nor end with a dot).

#. For hierarchies with virtual inheritance, document base virtual interfaces rather than overriding implementations.

#. Documented functions cannot use `GEOS_UNUSED_ARG()` in their declarations.

#. For empty virtual base implementations that use `GEOS_UNUSED_ARG(x)` to remove compiler warnings, use one of two options:

   * move empty definition away (e.g. out of class body) and keep `GEOS_UNUSED_ARG(x)` in definition only;
   * put `GEOS_UNUSED_VAR(x)` into the inline empty body.

#. For large classes, logically group functions using member groups via `///@{` and `///@}` and give them group names
   and descriptions (if needed) via a `@name` comment block. Typical groups may include:

   * constructors/destructor/assignment operators;
   * getter/setter type functions;
   * overridable virtual functions;
   * any other logically coherent groups (functions related to the same aspect of class behavior).

#. In-header implementation details (e.g. template helpers) often shouldn't appear in user documentation.
   Wrap these into `internal` namespace.

#. Use `/// @cond DO_NOT_DOCUMENT` and `/// @endcond` tags to denote a section of public API that should not be
   documented for some reason. This should be used rarely and selectively. An example is in-class helper structs
   that must be public but that user should not refer to explicitly.

Example
====================================

   .. code-block:: c++

      /// This is a documented macro
      #define USEFUL_MACRO

      /**
       * @brief Short description.
       * @tparam    T type of input value
       * @param[in] x input value explanation
       * @return      return value explanation
       *
       * Detailed description goes here.
       *
       * @note A note warning users of something unexpected.
       */
      template<typename T>
      int Foo( T const & x );

      /**
      * @brief Class for showing Doxygen.
      * @tparam T type of value the class operates on
      *
      * This class does nothing useful except show how to use Doxygen.
      */
      template<typename T>
      class Bar
      {
      public:

        /// A documented member type alias.
        using size_type = typename std::vector<T>::size_type;

        /**
         * @name Constructors/destructors.
         */
        ///@{

        /**
         * @brief A documented constructor.
         * @param value to initialize the object
         */
        explicit Bar( T t );

        /**
         * @brief A deleted, but still documented copy constructor.
         * @param an optionally documented parameter
         */
        Bar( Bar const & source ) = delete;

        /**
         * @brief A defaulted, but still documented move constructor.
         * @param an optionally documented parameter
         */
        Bar( Bar const & source ) = default;

        /**
         * @brief A documented desctructor.
         * virtual ~Bar() = default;
         */

        ///@}

        /**
         * @name Getters for stored value.
         */
        ///@{

        /**
         * @brief A documented public member function.
         * @return a reference to contained value
         */
        T & getValue();

        /**
         * @copydoc getValue()
         */
        T const & getValue() const;

        ///@}

      protected:

        /**
         * @brief A documented protected pure virtual function.
         * @param[in]  x the input value
         * @param[out] y the output value
         *
         * Some detailed explanation for users and implementers.
         */
        virtual void doSomethingOverridable( int const x, T & y ) = 0;

        /// @cond DO_NOT_DOCUMENT
        // Some stuff we don't want showing up in Doxygen
        struct BarHelper
        {};
        /// @endcond

      private:

        /// An optionally documented (not enforced) private member.
        T m_value;

      };

Current Doxygen
====================================

`Link to Doxygen Class directory <../../../../doxygen_output/html/classes.html>`__

Direct links to some useful class documentation:

`Group API <../../../../doxygen_output/html/classgeos_1_1data_repository_1_1_group.html>`_

`Wrapper API <../../../../doxygen_output/html/classgeos_1_1data_repository_1_1_wrapper.html>`_

`ObjectManagerBase API <../../../../doxygen_output/html/classgeos_1_1_object_manager_base.html>`_

`MeshLevel API <../../../../doxygen_output/html/classgeos_1_1_mesh_level.html>`_

`NodeManager API <../../../../doxygen_output/html/classgeos_1_1_node_manager.html>`_

`FaceManager API <../../../../doxygen_output/html/classgeos_1_1_face_manager.html>`_