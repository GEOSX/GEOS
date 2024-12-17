###############################################################################
Kernel interface
###############################################################################

Finite Element Method Kernel Interface
======================================

The finite element method kernel interface (FEMKI) specifies an API for the
launching of computational kernels for solving physics discretized using the
finite element method.
Using this approach, a set of generic element looping pattens and kernel
launching functions may be implemented, and reused by various physics solvers
that contain kernels conforming to the FEMKI.

There are several main components of the FEMKI:

#. A collection of element looping functions that provide various looping
   patterns, and call the ``launch`` function.

#. The kernel interface, which is specified by the
   `finiteElement::KernelBase <../../../doxygen_output/html/classgeos_1_1finite_element_1_1_implicit_kernel_base.html>`_ class.
   Each physics solver will define a class that contains its kernels functions,
   most likely deriving, or conforming to the API specified by the `KernelBase`
   class. Also part of this class will typically be a nested ``StackVariables``
   class that defines a collection of stack variables for use in the various
   kernel interface functions.

#. A ``launch`` function, which launches the kernel, and calls the kernel
   interface functions conforming to the interface defined by ``KernelBase``.
   This function is actually a member function of the ``Kernel`` class, so it
   may be overridden by a specific physics kernel, allowing complete
   customization of the interface, while maintaining the usage of the
   looping patterns.

A Generic Element Looping Pattern
---------------------------------
One example of a looping pattern is the
`regionBasedKernelApplication <../../../doxygen_output/html/_kernel_base_8hpp.html#file_a22560f1fcca889307fdabb5fa7422c0d>`_
function.

The contents of the looping function are displayed here:

.. literalinclude:: /coreComponents/finiteElement/kernelInterface/KernelBase.hpp
   :language: c++
   :start-after: //START_regionBasedKernelApplication
   :end-before: //END_regionBasedKernelApplication

This pattern may be used with any kernel class that either:

#. Conforms to the ``KernelBase`` interface by defining each of the kernel
   functions in ``KernelBase``.

#. Defines its own ``kernelLaunch`` function that conforms the the signature
   of ``KernelBase::kernelLaunch``.
   This option essentially allows for a custom kernel that does not conform to
   the interface defined by ``KernelBase`` and ``KernelBase::kernelLaunch``.

The KernelBase::kernelLaunch Interface
--------------------------------------
The ``kernelLaunch`` function is a member of the kernel class itself.
As mentioned above, a physics implementation may use the existing ``KernelBase``
interface, or define its own.
The ``KernelBase::kernelLaunch`` function defines a launching policy, and an
internal looping pattern over the quadrautre points, and calls the functions
defined by the ``KernelBase`` as shown here:

.. literalinclude:: /coreComponents/finiteElement/kernelInterface/KernelBase.hpp
   :language: c++
   :start-after: //START_kernelLauncher
   :end-before: //END_kernelLauncher

Each of the ``KernelBase`` functions called in the ``KernelBase::kernelLaunch``
function are intended to provide a certain amount of modularity and flexibility
for the physics implementations.
The general purpose of each function is described by the function name, but may
be further descibed by the function documentation found
`here <../../../doxygen_output/html/classgeos_1_1finite_element_1_1_kernel_base.html>`_.
