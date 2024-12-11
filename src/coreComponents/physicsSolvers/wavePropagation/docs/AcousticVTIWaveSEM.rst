
#####################################
Acoustic VTI Wave Propagation Solvers
#####################################

Available Solvers
=================


Governing equations
-------------------

The anisotropy is said to be TI when there exists an axis of symmetry that is normal to a plane of isotropy. A Vertical Transverse Isotropy (VTI) is a particular case of TI where the axis of symmetry correspond to the :math:`z`-axis or depth-axis. Numerous set of equations exist. Below are listed the available solvers in GEOS.

.. _sec-vti-fletcher:

``AcousticVTIFletcherWavePropagationSEM``
+++++++++++++++++++++++++++++++++++++++++

This solver relies on the set of equations proposed by Fletcher and his collaborators :footcite:`FlecherDuFowler2009`. This solver has a tunable parameter :math:`\vtif` (or equivalently :math:`\sigma`, see :ref:`paragraph on the quantities <sec-vti-quantities>`) to manage the spurious artefact s-wave.

.. math::
  :label: eqVTIFletcher

  \left\{
    \begin{aligned}
      \frac{1}{\rho v_p^2} \frac{\partial^2 p}{\partial t^2} &= (1+2\varepsilon)  \vtiAxy\nabla\cdot\left(\frac{1}{\rho} \vtiAxy\nabla p\right) + \vtiAz\nabla\frac{1}{\rho}\vtiAz\nabla q - (f_{\mathrm{vti}}-1)\vtiAz\nabla\cdot\left(\frac{1}{\rho}\vtiAz\nabla (p-q)\right) + f,\\
      \frac{1}{\rho v_p^2} \frac{\partial^2 q}{\partial t^2} &= (1+2\delta)  \vtiAxy\nabla\cdot\left(\frac{1}{\rho} \vtiAxy\nabla p\right)      + \vtiAz\nabla\frac{1}{\rho}\vtiAz\nabla q + (f_{\mathrm{vti}}-1) \vtiAxy\nabla\cdot\left(\frac{1}{\rho} \vtiAxy\nabla  (p-q)\right) + f.
   \end{aligned}
  \right.


.. _sec-vti-zhang:

``AcousticVTIZhangWavePropagationSEM``
++++++++++++++++++++++++++++++++++++++

This second one is based on Zhang et. al set of equations in :footcite:`ZhangZhangZhang2011`. As it is self-adjoint, this solver is mainly use for full waveform inversion.

.. math::
  :label: eqVTIZhang

  \left\{
    \begin{aligned}
       \frac{1}{\rho v_p^2} \frac{\partial^2 p}{\partial t^2} &= (1+2\varepsilon)  \vtiAxy\nabla\cdot\left(\frac{1}{\rho} \vtiAxy\nabla p\right) + \sqrt{1+2\delta} \vtiAz\nabla\cdot\left(\frac{1}{\rho}\vtiAz\nabla q\right) + f,\\
       \frac{1}{\rho v_p^2} \frac{\partial^2 q}{\partial t^2} &= \sqrt{1+2\delta}  \vtiAxy\nabla\cdot\left(\frac{1}{\rho} \vtiAxy\nabla p\right)      + \vtiAz\nabla\left(\frac{1}{\rho}\vtiAz\nabla q\right) + f.
    \end{aligned}
  \right.


.. _sec-vti-quantities:

Operators and quantities
------------------------

* :math:`p`: unknown pressure wave
* :math:`q`: auxiliary wave
* :math:`v_p [\textrm{m}.\textrm{s}^{-1}]`: pressure wave celerity in the medium.
* :math:`\rho [\textrm{kg}.\textrm{m}^{-3}]`: density of the medium.
* :math:`\varepsilon` and :math:`\delta`: anisotropy Thomsen parameters. Default value is 0 for both.
* :math:`\sigma_{\mathrm{vti}} = \frac{v_p^2}{v_s^2}(\varepsilon - \delta)`, with :math:`v_s` the share wave velocity, and :math:`f_{\mathrm{vti}} = 1 - \frac{\varepsilon - \delta}{\sigma_{\mathrm{vti}}}`: tunable parameters in Fletcher's equation to manage artefact spurious share wave. Even though :math:`\sigma_{\mathrm{vti}}` is the input in GEOS code, the equations are here written in term of :math:`f_{\mathrm{vti}}`. :math:`\sigma_{\mathrm{vti}}` is set by default to 0.75.
* :math:`f`: source term (Rickers)
* :math:`\nabla = [\partial x, \partial y, \partial z]^T`: nabla operator
* :math:`X\cdot Y = X^TY`: scalar product 
* :math:`\vtiAxy` and :math:`\vtiAz`: truncation matrices such that :math:`\vtiAxy\nabla = [\partial_x,\partial_y,0]^T` and :math:`\vtiAz\nabla=[0,0,\partial_z]^T`:

.. math::
  :label: eq-vti-AxyAz

  \vtiAxy = \begin{pmatrix}
  1&0&0\\
  0&1&0\\
  0&0&0
  \end{pmatrix}\quad\text{and}\quad
  \vtiAz = \begin{pmatrix}
  0&0&0\\
  0&0&0\\
  0&0&1
  \end{pmatrix}.

Restrictions
------------

- The computational domain must be a rectangular cuboid. Its boundary :math:`\Gamma` is divided into the lateral surfacess :math:`\Gammaxy` and the top/bottom surfaces :math:`\Gammaz`. The reason is to simply the boundary terms that arise in the weak formulation, as explained :ref:`in a next section <sec-vti-rectangularcuboid>`.
- The anisotropic parameters are (currently) constant per element.
- ``AcousticVTIFletcherWavePropagationSEM`` needs the following stability relation to be satisfied

  .. math::
    \varepsilon - \vtif^2 - \vtif\delta +\vtif + (1-\vtif)\sqrt{\vtif(\vtif+2\delta)} \geq 0.

- For numerical stability reasons, the ``AcousticVTIZhangWavePropagationSEM`` solver needs the following conditions to be satisfied:

  1. The :math:`\delta` parameter to be smooth in the domain (no sharp contrast). This must be achieved in the model by the user.
  2. :math:`\varepsilon \geq \delta` everywhere. If GEOS encounters the relation :math:`\delta > \varepsilon`, it will automatically set :math:`\delta = \varepsilon`.



Damping methods
---------------

.. _sec-vti-abc:

Absorbing Boundary Condition (ABC - default)
++++++++++++++++++++++++++++++++++++++++++++

The following ABC is set for both equations on the borders :math:`\Gamma_{xy}` and :math:`\Gamma_{z}` (See equation (2.66) in L. Boillot PhD thesis :footcite:`Boillot2014` )):

.. math::
  :label: eq-vti-alpha

  \begin{cases}
    \displaystyle \dn p = -\frac{1}{\alpha} \partial_t p,\\
    \displaystyle \dn q = -\frac{1}{\alpha} \partial_t q,
  \end{cases}\qquad
  \alpha = 
  \begin{cases}
    \displaystyle v_p\sqrt{1+2\varepsilon} &\Gamma_{xy},\\
    \displaystyle v_p &\Gamma_{z}.
  \end{cases}


.. _sec-vti-plm:

Perfectly Matched Layer (PML - not supported yet)
+++++++++++++++++++++++++++++++++++++++++++++++++

Perfectly Matched Layers (PML) are known to be unstable in a VTI media. They could however be used in conjunction with a taper to the model where the anisotropic parameters are vanishing (in a smooth way). The isotropic PML could then be applied as the medium is now isotropic at the border.

Currently, this option is not supported by GEOS though.

Parameters and Fields
---------------------

Supported Parameters
++++++++++++++++++++

The following attributes are supported, for the Fletcher's equations:

.. include:: /docs/sphinx/datastructure/AcousticVTIFletcherSEM.rst

And for the Zhang equations:

.. include:: /docs/sphinx/datastructure/AcousticVTIZhangSEM.rst


.. _sec-vti-fields:

Additional ``Fields``
+++++++++++++++++++++

The solvers use the exact same parameters as the isotropic one. The anisotropic parameters are implemented as GEOS ``Fields``. Another difference with isotropic solver is that the lateral and top/bottom surfaces must be properly defined. This is currently done using GEOS ``Fields``. 

.. list-table:: VTI Fields
   :header-rows: 1

   * - Name
     - Manager
     - Default
     - Description
   * - ``acousticEpsilon``
     - Cell
     - 0
     - :math:`\varepsilon` anisotropy Thomsen parameter
   * - ``acousticDelta``
     - Cell
     - 0
     - :math:`\delta` anisotropy Thomsen parameter
   * - ``acousticSigma``
     - Cell
     - 0.75
     - :math:`\sigma` anisotropy parameter (Fletcher's equation)
   * - ``acousticFreeSurface``
     - Face
     - 0
     - Set to 1 for nodes on the free surface (on top generally)
   * - ``acousticLateralSurface``
     - Face
     - 0
     - Set to 1 for nodes on the lateral surface :math:`\Gamma_{xy}`
   * - ``acousticBottomSurface``
     - Face
     - 0
     - Set to 1 for nodes on the bottom surface :math:`\Gamma_{z}`


.. _sec-vti-math:

Mathematical analysis
=====================

Explanations are here made for ``AcousticVTIZhangWavePropagationSEM`` with equations :eq:`eqVTIZhang`. The receipt is the same for ``AcousticVTIFletcherWavePropagationSEM`` and the results are summarized in the section :ref:`sec-vti-wf-final`.

Assumptions
-----------

* The anisotropic parameters are assumed to be constant per element. They will hence not be differentiated after the integration by part.
* The domain :math:`\Omega` is assumed to be a rectangular cuboid with boundary :math:`\Gamma = \Gamma_{xy} \bigcup \Gamma_{z}` where :math:`\Gamma_{xy}` is the lateral surface and :math:`\Gamma_{z}` represent the top and bottom surfaces.


General domain
--------------

After computation and gathering the unknown, the weak formulation reads as

.. math::
  \left\{
      \begin{aligned}
        &\text{Find } p,q\in C^2([0, +\infty])\times H^1(\Omega)  \text{ such that, }\forall p',q'\in H^1(\Omega)\times C^2([0, +\infty]),\\
        &\begin{multlined}[t]
           \int_{\Omega} \frac{1}{\rho v_p} \frac{\partial^2 p}{\partial t^2} p' \diff \mathbf{x} =
          - \int_{\Omega} \frac{(1+2\varepsilon)}{\rho} \vtiAxy\nabla p \cdot \vtiAxy\nabla p'\diff \mathbf{x}
          - \int_{\Omega} \frac{\sqrt{1+2\delta}}{\rho}\vtiAz\nabla q\cdot \vtiAz\nabla p'\diff \mathbf{x}\\
          + \int_{\Gamma} \frac{(1+2\varepsilon)}{\rho} (\vtiAxy\nabla p \cdot \vtiAxy \mathbf{n}) p'\diff s 
          + \int_{\Gamma} \frac{1}{\rho}\left(\vtiAz\nabla q\cdot \vtiAz\mathbf{n} \right)p'\diff s
          + \int_{\Omega} f p'\diff \mathbf{x},
        \end{multlined}\\
      &\begin{multlined}[t]
        \int_{\Omega} \frac{1}{\rho v_p} \frac{\partial^2 q}{\partial t^2} q' \diff \mathbf{x} = 
        - \int_{\Omega}\frac{\sqrt{1+2\delta}}{\rho} \vtiAxy\nabla p\cdot \vtiAxy\nabla q' \diff \mathbf{x} 
        - \int_{\Omega}\frac{1}{\rho}\vtiAz\nabla q\cdot \vtiAz\nabla q' \diff \mathbf{x}  \\
        + \int_{\Gamma}\frac{\sqrt{1+2\delta}}{\rho} (\vtiAxy\nabla p \cdot \vtiAxy\mathbf{n})  q' \diff s 
        + \int_{\Gamma}\frac{1}{\rho}\left(\vtiAz\nabla q\cdot \vtiAz\mathbf{n}\right) q' \diff s
        + \int_{\Omega}f  q'\diff \mathbf{x},
      \end{multlined}
    \end{aligned}
    \right.

.. _sec-vti-rectangularcuboid:

Rectangular Cuboid Domain (default)
-----------------------------------

This particular shape is handy for the boundary quantities. Indeed, on :math:`\Gammaz`, the normal vector :math:`\mathbf{n}` is equal to :math:`[0,0\pm 1]` and hence :math:`\vtiAz\mathbf{n} = \mathbf{n}`. The Neumann condition is recover as :math:`\vtiAz\nabla q\cdot \vtiAz\mathbf{n}`. The same idea can be applied on :math:`\Gammaxy` and the above set of equations can be slightly rewritten as

.. math::
  \left\{
      \begin{aligned}
        &\text{Find } p,q\in C^2([0, +\infty])\times H^1(\Omega)  \text{ such that, }\forall p',q'\in H^1(\Omega)\times C^2([0, +\infty]),\\
        &\begin{multlined}[t]
           \int_{\Omega} \frac{1}{\rho v_p} \frac{\partial^2 p}{\partial t^2} p' \diff \mathbf{x} =
          - \int_{\Omega} \frac{(1+2\varepsilon)}{\rho} \vtiAxy\nabla p \cdot \vtiAxy\nabla p'\diff \mathbf{x}
          - \int_{\Omega} \frac{\sqrt{1+2\delta}}{\rho}\vtiAz\nabla q\cdot \vtiAz\nabla p'\diff \mathbf{x}\\
          + \int_{\Gammaxy} \frac{(1+2\varepsilon)}{\rho} \dn(p) p'\diff s 
          + \int_{\Gammaz} \frac{1}{\rho} \dn(q)p'\diff s
          + \int_{\Omega} f p'\diff \mathbf{x},
        \end{multlined}\\
      &\begin{multlined}[t]
        \int_{\Omega} \frac{1}{\rho v_p} \frac{\partial^2 q}{\partial t^2} q' \diff \mathbf{x} = 
        - \int_{\Omega}\frac{\sqrt{1+2\delta}}{\rho} \vtiAxy\nabla p\cdot \vtiAxy\nabla q' \diff \mathbf{x} 
        - \int_{\Omega}\frac{1}{\rho}\vtiAz\nabla q\cdot \vtiAz\nabla q' \diff \mathbf{x}  \\
        + \int_{\Gammaxy}\frac{\sqrt{1+2\delta}}{\rho} \dn(p)  q' \diff s 
        + \int_{\Gammaz}\frac{1}{\rho}\dn(q) q' \diff s
        + \int_{\Omega}f  q'\diff \mathbf{x},
      \end{multlined}
    \end{aligned}
    \right.

GEOS assume that the shape of the computational domain is a rectangular cuboid.

Absorbing Boundary Conditions (ABC)
-----------------------------------

By default, the following ABC are set for both equations on the borders :math:`\Gamma_{xy}` and :math:`\Gamma_{z}` (See equation (2.66) in L. Boillot PhD thesis :footcite:`Boillot2014`):

.. math::
  :label: vti-alpha

  \begin{cases}
    \displaystyle \dn p = -\frac{1}{\alpha} \partial_t p,\\
    \displaystyle \dn q = -\frac{1}{\alpha} \partial_t q,
  \end{cases}\qquad
  \alpha = 
  \begin{cases}
    \displaystyle v_p\sqrt{1+2\varepsilon} &\Gamma_{xy},\\
    \displaystyle v_p &\Gamma_{z}.
  \end{cases}


Perfectly Matched Layer (PML)
-----------------------------

They are known to be unstable for VTI media. They could however be used in conjunction with a taper to the model where the anisotropic parameters are vanishing (in a smooth way). The PML can then be applied as the medium is now isotropic (on the border).

Currently, this option is not supported by GEOS.

.. _sec-vti-initial:

Initial condition
-----------------

Both waves and their derivatives are set to an initial value (default = 0):

.. math::

  \begin{cases}
    \displaystyle p(\vtixx,0) = p_0(\vtixx); \frac{\partial p}{\partial t}(\vtixx,0) = p_1(\vtixx), & \text{in }\Omega,\\
    \displaystyle q(\vtixx,0) = q_0(\vtixx); \frac{\partial q}{\partial t}(\vtixx,0) = q_1(\vtixx), & \text{in }\Omega.
  \end{cases}


Space discretization
--------------------

The unknown :math:`p` and :math:`q` are discretized using spectral element method or order :math:`r` leading to respectively the unknown vectors :math:`\vtipb` and :math:`\vtiqb` of :math:`\vtiVhr`. 
The following matrices are introduced where :math:`\Phi_I` and :math:`\Phi_J` refer to the basis functions associated to the :math:`I^{\textrm{th}}` and :math:`J^{\textrm{th}}` degree of freedom respectively. First, the mass and damping (or mass on the boundary) matrices

.. math::

  \left\{
  \begin{aligned}
    \vtiMass(\beta) &= \left(\vtiMass_{I,J}(\beta)\right)_{I,J},& \vtiMass_{I,J}(\beta) & = \vtiint{\Omega}{\beta(\vtixx)\Phi_J(\vtixx)\Phi_I(\vtixx)}{\vtixx},\\
    \vtiDamp(\beta) &= \left(\vtiDamp_{I,J}(\beta)\right)_{I,J},& \vtiDamp_{I,J} (\beta)& = \vtiint{\Gamma}{\beta(s(\vtixx))\Phi_J(s(\vtixx))\Phi_I(s(\vtixx))}{s},\\
    \vtiDampxy(\beta) &= \left(\vtiDamp^{xy}_{I,J}(\beta)\right)_{I,J},& \vtiDamp^{xy}_{I,J}(\beta) & = \vtiint{\Gammaxy}{\beta(s(\vtixx))\Phi_J(s(\vtixx))\Phi_I(s(\vtixx))}{s},\\
    \vtiDampz(\beta) &= \left(\vtiDamp^{z}_{I,J}(\beta)\right)_{I,J},& \vtiDamp^{z}_{I,J} (\beta)& = \vtiint{\Gammaz}{\beta(s(\vtixx))\Phi_J(s(\vtixx))\Phi_I(s(\vtixx))}{s}.
  \end{aligned}
  \right.

Second, the stiffness and generalized stiffness matrices are defined by

.. math::
  :label: eq-vti-stiff

  \left\{
  \begin{aligned}
    \vtiStiff(\beta) &=  \left(\vtiStiff_{I,J}(\beta)\right)_{I,J},& \vtiStiff_{I,J} (\beta)& = \vtiint{\Omega}{\beta(\vtixx)\nabla \Phi_J(\vtixx)\cdot\nabla\Phi_I(\vtixx)}{\vtixx},\\
    \vtiStiffxy(\beta) &= \left(\vtiStiff^{xy}_{I,J}(\beta)\right)_{I,J},& \vtiStiff^{xy}_{I,J} (\beta)& = \vtiint{\Omega}{\beta(\vtixx)\vtiAxy\Phi_J(\vtixx)\cdot\vtiAxy\nabla\Phi_I(\vtixx)}{\vtixx},\\
    \vtiStiffz(\beta) &= \left(\vtiStiff^z_{I,J}(\beta)\right)_{I,J},& \vtiStiff^z_{I,J} (\beta)& = \vtiint{\Omega}{\beta(\vtixx)\vtiAz\nabla\Phi_J(\vtixx)\cdot\vtiAz\nabla \Phi_I(\vtixx)}{\vtixx}.
  \end{aligned}
  \right.



The discretized weak formulation is then given by

.. 
  .. math::
    \left\{
      \begin{aligned}
        &\text{Find } \vtipb,\vtiqb\in \vtiVhr\times C^2([0, +\infty])  \text{ such that, }\forall \vtipb',\vtiqb'\in \vtiVhr\times C^2([0, +\infty]),\\
        &\begin{multlined}[t]
          \frac{1}{\rho v_p^2}\vtiMass \frac{\partial^2\vtipb}{\partial t^2} =
          -\vtiStiffxy(1+2\varepsilon) \vtipb 
          - \vtiStiffz \vtiqb
          + \vtiStiffz(\vtif-1)(\vtipb-\vtiqb) \\
          - \vtiDampxy(\alpha(1+2\varepsilon)) \frac{\partial \vtipb}{\partial t} 
          - \vtiDampz(\alpha)\frac{\partial \vtiqb}{\partial t}
          + \vtiDampz(\alpha(\vtif-1))\frac{\partial (\vtipb-\vtiqb)}{\partial t}
          + \vtiMass\vtifbq,
        \end{multlined}\\
      &\begin{multlined}[t]
        \frac{1}{\rho v_p^2}\vtiMass \frac{\partial^2\vtiqb}{\partial t^2} = 
        -\vtiStiffxy(1+2\delta) \vtipb 
        - \vtiStiffz \vtiqb
        - \vtiStiffxy(\vtif-1) (\vtipb-\vtiqb)\\
        -\vtiDampxy(\alpha(1+2\delta)) \frac{\partial\vtipb}{\partial t} 
        - \vtiDampz(\alpha)\frac{\partial\vtiqb}{\partial t}
        - \vtiDampxy(\alpha(\vtif-1)) \frac{\partial(\vtipb-\vtiqb)}{\partial t}
        + \vtiMass\vtifbq.
      \end{multlined}
    \end{aligned}
    \right.


.. math::

  \left\{
    \begin{aligned}
      &\begin{multlined}
        \vtiMass\left(\frac{1}{\rho v_p^2}\right)\frac{\partial^2 \vtipb}{\partial t^2}
        + \vtiStiffxy\left(\frac{1+2\varepsilon}{\rho}\right) \vtipb
        + \vtiStiffz\left(\frac{\sqrt{1+2\delta}}{\rho}\right) \vtiqb \\
        + \vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)\frac{\partial \vtipb}{\partial t}
        + \vtiDampz\left( \frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \frac{\partial \vtiqb}{\partial t}
        = \vtifbp,
      \end{multlined}\\
      &\begin{multlined}
        \vtiMass\left(\frac{1}{\rho v_p^2}\right) \frac{\partial^2 \vtiqb }{\partial t^2}
        + \vtiStiffxy\left(\frac{\sqrt{1+2\delta}}{\rho}\right) \vtipb
        + \vtiStiffz \left(\frac{1}{\rho}\right) \vtiqb\\
        + \vtiDampxy\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right)\frac{\partial \vtipb}{\partial t}
        + \vtiDampz \left(\frac{\alpha                }{\rho}\right)\frac{\partial \vtiqb}{\partial t}
        = \vtifbq.
      \end{multlined}
    \end{aligned}
  \right.


Time discretization
-------------------

The leapfrog scheme with :math:`\vtidt` as a time step leads to the following approximation to compute time step :math:`n+1` from :math:`n` and :math:`n-1`:

.. math::
  \left\{
    \begin{aligned}
      \frac{\partial^2\vtipb}{\partial t^2} &\approx  \frac{\vtipb^{n+1} - 2\vtipb^{n} +\vtipb^{n-1}}{\vtidt^2},\\
      \frac{\partial \vtipb}{\partial t} &\approx \frac{\vtipb^{n+1}- \vtipb^{n-1}}{2 \vtidt},
    \end{aligned}
  \right.
  \quad\text{ and }\quad
  \left\{
    \begin{aligned}
      \frac{\partial^2\vtiqb}{\partial t^2} &\approx  \frac{\vtiqb^{n+1} - 2\vtiqb^{n} +\vtiqb^{n-1}}{\vtidt^2},\\
      \frac{\partial\vtiqb}{\partial t} &\approx \frac{\vtiqb^{n+1}- \vtiqb^{n-1}}{2 \vtidt}.
    \end{aligned}
  \right.

.. _sec-vti-wf-final:

Weak Formulations (Final Form)
------------------------------


The weak formulation for both solvers finally read as followm where the ABC parameter :math:`\alpha` given by equation :eq:`vti-alpha`.


Weak formulation for ``AcousticVTIFletcherWavePropagationSEM``
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. math::

  \left\{
    \begin{aligned}
      &\text{Find } \vtipb,\vtiqb\in \mathbb{R}^{d_r}  \text{ such that,}\\
      &\begin{multlined}[t]
        \left[\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)
        +\frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)
        - \frac{1}{2\vtidt}\vtiDampz\left(\frac{\alpha(\vtif-1)}{\rho}\right)
        \right]\vtipb^{n+1} 
        + \frac{1}{2\vtidt} \vtiDampz\left(\frac{\alpha \vtif}{\rho}\right) \vtiqb^{n+1}=  
         \frac{2}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)  \vtipb^{n} 
        -\frac{1}{ \vtidt^2}   \vtiMass\left(\frac{1}{\rho v_p^2}\right)\vtipb^{n-1} \\
        -\vtiStiffxy(1+2\varepsilon) \vtipb^n 
        + \vtiStiffz(\vtif-1)\vtipb^n
        -\vtiStiffz(\vtif)\vtiqb^n 
        + \left[\frac{1}{2 \vtidt}\vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)
        - \frac{1}{2 \vtidt}\vtiDampz\left(\frac{\alpha(\vtif-1))}{\rho}\right)\right]\vtipb^{n-1}
        +\frac{1}{2 \vtidt}\vtiDampz\left(\frac{\alpha \vtif}{\rho}\right)\vtiqb^{n-1}
        + \vtiMass\vtifbp^n,
      \end{multlined}\\
      &\begin{multlined}[t]
        \left[\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right) 
        + \frac{1}{2\vtidt}\vtiDampz \left(\frac{\alpha}{\rho}\right)
        - \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha(\vtif-1)}{\rho}\right)\right]\vtiqb^{n+1}
        + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha (\vtif + 2\delta)}{\rho}\right) \vtipb^{n+1}
        =  \frac{1}{2\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2 }\right) \vtiqb^{n} -\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2 }\right)\vtiqb^{n-1} \\
        -\vtiStiffxy(2\delta+\vtif) \vtipb^n 
        + \vtiStiffxy(\vtif-1) \vtiqb^n
        - \vtiStiffz \vtiqb^n
        +\frac{1}{2 \vtidt}\vtiDampxy\left(\frac{\alpha(2\delta + \vtif)}{\rho}\right) \vtipb^{n-1}
        +\frac{1}{2 \vtidt}\left[- \vtiDampxy\left(\frac{\alpha(\vtif-1)}{\rho}\right)
        +\frac{1}{2 \vtidt}\vtiDampz\left(\frac{\alpha}{\rho}\right) \right]\vtiqb^{n-1}
        + \vtiMass\vtifbq^n.
      \end{multlined}
    \end{aligned}
  \right.



Weak formulation for ``AcousticVTIZhangWavePropagationSEM``
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


.. math::
  :label: eq-vti-zhang-wf
  
  \left\{\begin{aligned}
    &\begin{multlined}
    \left[\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right) 
    + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)\right]\vtipb^{n+1} 
    + \frac{1}{2\vtidt}\vtiDampz\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \vtiqb^{n+1} =
     \left[\frac{2}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right) 
    - \vtiStiffxy\left(\frac{1+2\varepsilon  }{\rho}\right)\right] \vtipb^n\\
    - \vtiStiffz \left(\frac{\sqrt{1+2\delta}}{\rho}\right) \vtiqb^n 
    + \left[-\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)
    + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)\right] \vtipb^{n-1}
    + \frac{1}{2\vtidt}\vtiDampz \left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \vtiqb^{n-1}
    + \vtifbp,
    \end{multlined}\\
    &\begin{multlined}
  \left[\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right) 
    + \frac{1}{2\vtidt}\vtiDampz \left(\frac{\alpha                }{\rho}\right)\right]\vtiqb^{n+1}
    + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \vtipb^{n+1} =
    \left[\frac{2}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)
    - \vtiStiffz \left(\frac{1               }{\rho}\right)\right]\vtiqb^n \\
    - \vtiStiffxy\left(\frac{\sqrt{1+2\delta}}{\rho}\right)\vtipb^n
    + \left[-\frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)
    +\frac{1}{2\vtidt}\vtiDampz  \left(\frac{\alpha                }{\rho}\right) \right] \vtiqb^{n-1}
    + \frac{1}{2\vtidt}\vtiDampxy\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right)\vtipb^{n-1}
    + \vtifbq.
    \end{multlined}
    \end{aligned}\right.


Implementation
==============

.. _sec-vti-gen-stiff:

Generalized Stiffness Matrices
------------------------------

General case
++++++++++++


In this section is explained how are computed the terms :math:`\vtiStiffxy` and :math:`\vtiStiffz` in equation :eq:` eq-vti-stiff` and actually how, in GEOS, are computed every general stiffness matrix :math:`\vtiStiff^{\mathbf{A}}` for a matrix :math:`\mathbf{A}`:

.. math::
  :label: eq-vti-stiffA

  \vtiStiff^{\mathbf{A}}_{\vtiib,\vtijb} = \vtiint{\Omega}{\left(\mathbf{A}\nabla \Phi_{\vtijb}\right)\cdot\left(\mathbf{A}\nabla\Phi_{\vtiib}\right)}{\vtixx}.

For the sake of clarity, every Degree of Freedom (DoF) of the SEM of order :math:`r` are assumed to be numbered using a triplet :math:`\vtiib=\{i_1,i_2,i_3\} \in \mathbb{I}`, with

.. math::
  \mathbb{I} = [ 1, 2, \ldots, r+1]^3 = \vtienstq{ \vtiib=(i_1,i_2,i_3) \in \mathbb{N}^3}{1\leq i_1,i_2,i_3 \leq r+1}.

The basis function associated to the DoF number :math:`\vtiib` is here denoted by :math:`\Phi_{\vtiib}`. The mesh is composed by :math:`N_{\textrm{elem}}` hexahedra denoted :math:`K_p` for :math:`p=1,2,\ldots,N_{\textrm{elem}}`. 
The (classical) stiffness coefficients  :math:`\vtiStiff_{\vtiib, \vtijb}` is obtained through the addition of every elementary contribution:

.. math::
  \vtiStiff^p_{\vtiib, \vtijb} = \int_{K_p} \nabla\Phi_{\vtijb}(\vtixx) \cdot\nabla\Phi_{\vtiib}(\vtixx)\diff\vtixx.

Computing :math:`\vtiStiff^p_{\vtiib, \vtijb}` is done through a transformation, with jacobian matrix :math:`\vtiJacp`, from the reference element :math:`\hat{K}` to the current one:

.. math::
  \vtiStiff^p_{\vtiib, \vtijb} = \int_{\hat{K}} \left(\vtiBbp(\hat{\vtixx})\nabla\hat{\Phi}_{\vtijb}(\hat{\vtixx})\right) \cdot\nabla\hat{\Phi}_{\vtiib}(\hat{\vtixx})\diff\hat{\vtixx},
  \quad\text{with}\quad
  \vtiBbp(\hat{\vtixx}) = \abs{\vtijacp(\hat{\vtixx})}\vtiJacp^{-1}(\hat{\vtixx})\vtiJacp^{-T}(\hat{\vtixx}).

.. math::
  \vtiStiff^p_{\vtiib, \vtijb} = \int_{\hat{K}} \left(\vtiBbp(\hat{\vtixx})\nabla\hat{\Phi}_{\vtijb}(\hat{\vtixx})\right) \cdot\nabla\hat{\Phi}_{\vtiib}(\hat{\vtixx})\diff\hat{\vtixx},
  \quad\text{with}\quad
  \vtiBbp(\hat{\vtixx}) = \abs{\vtijacp(\hat{\vtixx})}\vtiJacp^{-1}(\hat{\vtixx})\vtiJacp^{-T}(\hat{\vtixx}).


Now, the same procedure is applied to :math:`\vtiStiff^{\mathbf{A}}` from equation :eq:`eq-vti-stiffA` and, to simplify, the quantities :math:`\vtixx` are hidden. The elementary contribution by element :math:`K_p` is defined by

.. math::
  \vtiStiff^{p,\mathbf{A}}_{\vtiib, \vtijb} = \int_{K_p} \left(\mathbf{A}\nabla\Phi_{\vtijb}\right) \cdot\left(\mathbf{A} \nabla\Phi_{\vtiib}\right)\diff\vtixx,

and after transformation from the reference element, we obtain

.. math::
  \vtiStiff^{p,\mathbf{A}}_{} = \int_{\hat{K}} \left(\vtiBbpA(\hat{\vtixx})\nabla\hat{\Phi}_{\vtijb}(\hat{\vtixx})\right) \cdot\nabla\hat{\Phi}_{\vtiib}(\hat{\vtixx})\diff\hat{\vtixx},

where the matrix  :math:`\vtiBbpA` is given by

.. math::
  \vtiBbpA(\hat{\vtixx}) = \abs{\vtijacp(\hat{\vtixx})}\vtiJacp^{-1}(\hat{\vtixx})\mathbf{A(\hat{\vtixx})}^T\mathbf{A(\hat{\vtixx})}\vtiJacp^{-T}(\hat{\vtixx}).

The only difference with the standard stiffness is thus the :math:`\vtiBbp` matrix.


Case :math:`\mathbf{A}=\vtiAxy` or :math:`\mathbf{A}=\vtiAz`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

These matrices, given by :eq:`eq-vti-AxyAz`, satisfy :math:`\mathbf{A}^T\mathbf{A}=\mathbf{A}`. Their associated matrices :math:`\vtiBbpxy` and :math:`\vtiBbpz` can be easily computed and are given by


.. math::

  \begin{multlined}
    \vtiBbpxy = \abs{\vtijacp}
    \left(\begin{matrix*}[l]
    \vtiJacp^{-1}[1,1]^2 + \vtiJacp^{-1}[1,2]^2  & \vtiJacp^{-1}[1,1]\vtiJacp^{-1}[2,1] +\vtiJacp^{-1}[1,2]\vtiJacp^{-1}[2,2] \\
    \vtiJacp^{-1}[2,1]\vtiJacp^{-1}[1,1] + \vtiJacp^{-1}[2,2]\vtiJacp^{-1}[1,2] & \vtiJacp^{-1}[2,2]^2 + \vtiJacp^{-1}[2,1]^2 \\
    \vtiJacp^{-1}[3,1]\vtiJacp^{-1}[1,1]\vtiJacp^{-1}[3,2]\vtiJacp^{-1}[1,2] & \vtiJacp^{-1}[3,1]\vtiJacp^{-1}[2,1]\vtiJacp^{-1}[3,2]\vtiJacp^{-1}[2,2] 
    \end{matrix*}\right.\\
    \left.\begin{matrix*}[l]
      \vtiJacp^{-1}[1,1]\vtiJacp^{-1}[3,1] + \vtiJacp^{-1}[1,2]\vtiJacp^{-1}[3,2]\\
      \vtiJacp^{-1}[2,1]\vtiJacp^{-1}[3,1]+\vtiJacp^{-1}[2,2]\vtiJacp^{-1}[3,2]\\
      \vtiJacp^{-1}[3,1]^2+ \vtiJacp^{-1}[3,2]^2
    \end{matrix*}\right),
  \end{multlined}

and 

.. math::
  \vtiBbpz = \abs{\vtijacp}
  \begin{pmatrix*}[l]
    \vtiJacp^{-1}[1,3]^2 & \vtiJacp^{-1}[1,3]\vtiJacp^{-1}[2,3] & \vtiJacp^{-1}[1,3]\vtiJacp^{-1}[3,3]\\
    \vtiJacp^{-1}[2,3]\vtiJacp^{-1}[1,3] & \vtiJacp^{-1}[2,3]^2 & \vtiJacp^{-1}[2,3]\vtiJacp^{-1}[3,3]\\
    \vtiJacp^{-1}[3,3]\vtiJacp^{-1}[1,3] & \vtiJacp^{-1}[3,3]\vtiJacp^{-1}[2,3] & \vtiJacp^{-1}[3,3]^2
  \end{pmatrix*}.

The computation of :math:`\vtiStiffxy` and :math:`\vtiStiffz` should then be about the same complexity as the computation of the stiffness matrix :math:`\vtiStiff`.

.. _sec-vti-boundary-nodes:

Interior vs Boundary Nodes
--------------------------

In SEM and contrary to classical FEM, the mass matrix is diagonal and so are the damping matrices too. For an interior degree of freedom, numbered :math:`\vtiib`, computing the next time step :math:`\vtipb^{n+1}_{\vtiib}` from :eq:`eq-vti-zhang-wf` reduces to a simple division

.. math::
  \left\{\begin{aligned}
    \frac{1}{\vtidt^2}\vtiMass_{\vtiib,\vtiib}\left(\frac{1}{\rho v_p^2}\right) \vtipb_{\vtiib}^{n+1} 
    &=
    a_{\vtiib},\\
    \frac{1}{\vtidt^2}\vtiMass\left(\frac{1}{\rho v_p^2}\right)  \vtiqb_{\vtiib}^{n+1}
     &= 
     b_{\vtiib}.
    \end{aligned}\right.

The quantities :math:`a_{\vtiib}` and :math:`b_{\vtiib}` are the :math:`\vtiib^{\textrm{th}}` component of the right hand sides of :eq:`eq-vti-zhang-wf`. For an interior node however, a coupling between :math:`\vtipb_{\vtiib}^{n+1}` and :math:`\vtiqb_{\vtiib}^{n+1}` appears:

.. math::
  \left\{\begin{aligned}
    \left[\frac{1}{\vtidt^2}\vtiMass_{\vtiib,\vtiib}\left(\frac{1}{\rho v_p^2}\right) 
    + \frac{1}{2\vtidt}\vtiDampxyii\left(\frac{\alpha(1+2\varepsilon)}{\rho}\right)\right]\vtipb_{\vtiib}^{n+1} 
    + \frac{1}{2\vtidt}\vtiDampzii\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \vtiqb_{\vtiib}^{n+1}
    &=a_{\vtiib}\\
      \left[\frac{1}{\vtidt^2}\vtiMass_{\vtiib,\vtiib}\left(\frac{1}{\rho v_p^2}\right) 
      + \frac{1}{2\vtidt}\vtiDampzii \left(\frac{\alpha                }{\rho}\right)\right]\vtiqb_{\vtiib}^{n+1}
      + \frac{1}{2\vtidt}\vtiDampxyii\left(\frac{\alpha\sqrt{1+2\delta}}{\rho}\right) \vtipb_{\vtiib}^{n+1}
     &= 
     b_{\vtiib}
    \end{aligned}\right.

This is simple :math:`2\times 2` system with unknown :math:`[\vtipb_{\vtiib}^{n+1},\vtiqb_{\vtiib}^{n+1}]^T` which is solved analytically in GEOS for each degree of freedom at the boundary.


.. footbibliography::
