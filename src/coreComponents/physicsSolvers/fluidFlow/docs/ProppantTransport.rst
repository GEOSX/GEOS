.. _ProppantTransport:

#####################################
Proppant Transport Solver
#####################################

Introduction
=========================

The `ProppantTransport` solver applies the finite volume method to solve the equations of proppant transport in hydraulic fractures. The behavior of proppant transport is described by a continuum formulation. Here we briefly outline the usage, governing equations and numerical implementation of the proppant transport model in GEOS.

Theory
=========================

The following mass balance and constitutive equations are solved inside fractures,

Proppant-fluid Slurry Flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
   \frac{\partial}{\partial t}(\rho_m) + \boldsymbol{\nabla} \cdot (\rho_m \boldsymbol{u_m}) = 0,

where the proppant-fluid mixture velocity :math:`\boldsymbol{u_m}` is approximated by the Darcy's law as,

.. math::
   \boldsymbol{u}_m = -\frac{K_f}{\mu_m}(\nabla p - \rho_m \boldsymbol{g}),

and :math:`p` is pressure, :math:`\rho_m` and :math:`\mu_m` are density and viscosity of the mixed fluid , respectively,  and :math:`\boldsymbol{g}` is the gravity vector. The fracture permeability :math:`K_f` is determined based on fracture aperture :math:`a` as

.. math::
   K_f =  \frac{a^2}{12}


Proppant Transport
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
   \frac{\partial}{\partial t}(c) + \boldsymbol{\nabla} \cdot (c \boldsymbol{u}_p) = 0,

in which :math:`c` and :math:`\boldsymbol{u}_p` represent the volume fraction and velocity of the proppant particles.


Multi-component Fluid Transport
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
   \frac{\partial}{\partial t} [ \rho_i \omega_i (1 - c) ] + \boldsymbol{\nabla} \cdot [ \rho_i \omega_i (1 - c) \boldsymbol{u}_f ] = 0.

Here :math:`\boldsymbol{u}_f` represents the carrying fluid velocity. :math:`\rho_i` and :math:`\omega_i` denote the density and concentration of `i-th` component in fluid, respectively. The fluid density :math:`\rho_f` can now be readily written as

.. math::
   \rho_f = \sum_{i=1}^{N_c} \rho_i \omega_i,

where :math:`N_c` is the number of components in fluid.
Similarly, the fluid viscosity :math:`\mu_f` can be calculated by the mass fraction weighted average of the component viscosities.

The density and velocity of the slurry fluid are further expressed as,

.. math::
   \rho_m = (1 - c) \rho_f + c \rho_p,

and

.. math::
   \rho_m \boldsymbol{u}_m = (1 - c) \rho_f \boldsymbol{u}_f + c \rho_p \boldsymbol{u}_p,

in which :math:`\rho_f` and :math:`\boldsymbol{u}_f` are the density and velocity of the carrying fluid, and :math:`\rho_p` is the density of the proppant particles.


Proppant Slip Velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The proppant particle and carrying fluid velocities are related by the slip velocity :math:`\boldsymbol{u}_{slip}`,

 .. math::
    \boldsymbol{u}_{slip} = \boldsymbol{u}_p - \boldsymbol{u}_f.

The slip velocity between the proppant and carrying fluid includes gravitational and collisional components, which take account of particle settling and collision effects, respectively.

The gravitational component of the slip velocity :math:`\boldsymbol{u}_{slipG}` is written as a form as

.. math::
    \boldsymbol{u}_{slipG} = F(c) \boldsymbol{u}_{settling},


where :math:`\boldsymbol{u}_{settling}` is the settling velocity for a single particle, :math:`d_p` is the particle diameter, and :math:`F(c)` is the correction factor to the particle settling velocity in order to account for hindered settling effects as a result of particle-particle interactions,

.. math::
    F(c) = e^{-\lambda_s c},

with the hindered settling coefficient :math:`\lambda_s` as an empirical constant set to 5.9 by default (Barree & Conway, 1995).

The settling velocity for a single particle, :math:`\boldsymbol{u}_{settling}` , is calculated based on the Stokes drag law by default,

.. math::
    \boldsymbol{u}_{settling} = ( \rho_p - \rho_f)  \frac{d{_p}^{2}}{18 \mu_f}\boldsymbol{g}.

Single-particle settling under intermediate Reynolds-number and turbulent flow conditions can also be described respectively by the Allen's equation (Barree & Conway, 1995),

.. math::
    \boldsymbol{u}_{settling} = 0.2 d_{p}^{1.18} \left [ \frac{g ( \rho_p - \rho_f)}{\rho_f} \right ]^{0.72} \left ( \frac{\rho_f}{\mu_f} \right )^{0.45} \boldsymbol{e},

and Newton's equation(Barree & Conway, 1995),

.. math::
    \boldsymbol{u}_{settling} = 1.74 d{_p}^{0.5}\left [ \frac{g ( \rho_p - \rho_f)}{\rho_f}\right]^{0.5} \boldsymbol{e}.


:math:`\boldsymbol{e}` is the unit gravity vector and :math:`d_p` is the particle diameter.

The collisional component of the slip velocity is modeled by defining :math:`\lambda`, the ratio of the particle velocity to the volume averaged mixture velocity as a function of the proppant concentration. From this the particle slip velocity in horizontal direction is related to the mixed fluid velocity by,

.. math::
    \boldsymbol{u}_{slipH} =  \frac{\lambda - 1}{1 - c} \boldsymbol{v}_{m}

with :math:`\boldsymbol{v}_{m}` denoting volume averaged mixture velocity.
We use a simple expression of :math:`\lambda` proposed by Barree & Conway (1995) to correct the particle slip velocity in horizontal direction,

.. math::
    \lambda=  \left[\alpha - |c - c_{slip} |^{\beta} \right]\,

where :math:`\alpha` and :math:`\beta` are empirical constants, :math:`c_{slip}` is the volume fraction exhibiting the greatest particle slip. By default the model parameters are set to the values given in (Barree & Conway, 1995): :math:`\alpha= 1.27`, :math:`c_{slip} =0.1` and :math:`\beta =  1.5`. This model can be extended to account for the transition to the particle pack as the proppant concentration approaches the jamming transition.

Proppant Bed Build-up and Load Transport
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to suspended particle flow the GEOS has the option to model proppant settling into an immobile bed at the bottom of the fracture. As the proppant cannot settle further down the proppant bed starts to form and develop at the element that is either at the bottom of the fracture or has an underlying element already filled with particles. Such an "inter-facial" element is divided into proppant flow and immobile bed regions based on the proppant-pack height.

Although proppant becomes immobile fluid can continue to flow through the settled proppant pack. The pack permeability `K` is defined based on the Kozeny-Carmen relationship:

.. math::
   K = \frac{(sd_p)^2}{180}\frac{\phi^{3}}{(1-\phi)^{2}}

and

.. math::
   \phi = 1 - c_{s}

where :math:`\phi` is the porosity of particle pack and :math:`c_{s}` is the saturation or maximum fraction for proppant packing, :math:`s` is the sphericity and :math:`d_p` is the particle diameter.


The growth of the settled pack in an "inter-facial" element is controlled by the interplay between proppant gravitational settling and shear-force induced lifting as (Hu et al., 2018),

.. math::
    \frac{d H}{d t} =  \frac{c u_{settling} F(c)}{c_{s}} - \frac{Q_{lift}}{A c_{s}},

where :math:`H`, :math:`t`, :math:`c_{s}`, :math:`Q_{lift}`, and :math:`A` represent the height of the proppant bed, time, saturation or maximum proppant concnetration in the proppant bed, proppant-bed load (wash-out) flux, and cross-sectional area, respectively.

The rate of proppant bed load transport (or wash out) due to shear force is calculated by the correlation proposed by Wiberg and Smith (1989) and McClure (2018),

.. math::
    Q_{lift} = a \left ( d{_p} \sqrt{\frac{g d{_p} ( \rho_p - \rho_f)}{\rho_f}} \right ) (9.64 N_{sh}^{0.166})(N_{sh} - N_{sh, c})^{1.5}.

:math:`a` is fracture aperture, and :math:`N_{sh}` is the Shields number measuring the relative importance of the shear force to the gravitational force on a particle of sediment (Miller et al., 1977; Biot & Medlin, 1985; McClure, 2018) as

.. math::
   N_{sh} = \frac{\tau}{d{_p} g ( \rho_p - \rho_f)},

and

.. math::
   \tau = 0.125 f \rho_f u_{m}^2

where :math:`\tau` is the shear stress acting on the top of the proppant bed and :math:`f` is the Darcy friction coefficient. :math:`N_{sh, c}` is the critical Shields number for the onset of bed load transport.


Proppant Bridging and Screenout
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Proppant bridging occurs when proppant particle size is close to or larger than fracture aperture. The aperture at which bridging occurs, :math:`h_{b}`, is defined simply by

.. math::
   h_{b} = \lambda_{b} d_p,

in which :math:`\lambda_{b}` is the bridging factor.

Slurry Fluid Viscosity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The viscosity of the bulk fluid, :math:`\mu_m`, is calculated as a function of proppant concentration as (Keck et al., 1992),

.. math::
     \mu_{m} =  \mu_{f}\left [1 + 1.25 \left ( \frac{c}{1-c/c_{s}} \right) \right ]^{2}.


Note that continued model development and improvement are underway and additional empirical correlations or functions will be added to support the above calculations.


Spatial Discretization
=======================

The above governing equations are discretized using a cell-centered two-point flux approximation (TPFA) finite volume method. We use an upwind scheme to approximate proppant and component transport across cell interfaces.


Solution Strategy
=======================

The discretized non-linear slurry flow and proppant/component transport equations at each time step are separately solved by the Newton-Raphson method. The coupling between them is achieved by a time-marching sequential (operator-splitting) solution approach.


Parameters
=========================

The solver is enabled by adding a ``<ProppantTransport>`` node
and a ``<SurfaceGenerator>`` node in the Solvers section.
Like any solver, time stepping is driven by events, see :ref:`EventManager`.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/ProppantTransport.rst

In particular:

* ``discretization`` must point to a Finite Volume flux approximation scheme defined in the Numerical Methods section of the input file (see :ref:`FiniteVolume`)
* ``proppantName`` must point to a particle fluid model defined in the Constitutive section of the input file (see :ref:`Constitutive`)
* ``fluidName`` must point to a slurry fluid model defined in the Constitutive section of the input file (see :ref:`Constitutive`)
* ``solidName`` must point to a solid mechanics model defined in the Constitutive section of the input file (see :ref:`Constitutive`)
* ``targetRegions`` attribute is currently not supported, the solver is always applied to all regions.

Primary solution field labels are ``proppantConcentration`` and
``pressure``.
Initial conditions must be prescribed on these field in every region, and boundary conditions
must be prescribed on these fields on cell or face sets of interest. For static (non-propagating) fracture problems, the fields ``ruptureState`` and
``elementAperture`` should be provided in the initial conditions.

In addition, the solver declares a scalar field named ``referencePorosity`` and a vector field
named ``permeability``, that contains principal values of the symmetric rank-2 permeability tensor
(tensor axis are assumed aligned with the global coordinate system).
These fields must be populated via :ref:`XML_FieldSpecification` section and ``permeability`` should
be supplied as the value of ``coefficientName`` attribute of the flux approximation scheme used.

Example
=========================

First, we specify the proppant transport solver itself and apply it to the fracture region:

.. literalinclude:: ../../../../../inputFiles/proppant/FlowProppantTransport2d_base.xml
  :language: xml
  :start-after: <!-- SPHINX_PROPPANT_TRANSPORT_SOLVER_BEGIN -->
  :end-before: <!-- SPHINX_PROPPANT_TRANSPORT_SOLVER_END -->

Then, we specify a compatible flow solver (currently a specialized ``SinglePhaseProppantFVM`` solver must be used):

.. literalinclude:: ../../../../../inputFiles/proppant/FlowProppantTransport2d_base.xml
  :language: xml
  :start-after: <!-- SPHINX_PROPPANT_FLOW_SOLVER_BEGIN -->
  :end-before: <!-- SPHINX_PROPPANT_FLOW_SOLVER_END -->

Finally, we couple them through a coupled solver that references the two above:

.. literalinclude:: ../../../../../inputFiles/proppant/FlowProppantTransport2d_base.xml
  :language: xml
  :start-after: <!-- SPHINX_PROPPANT_COUPLED_SOLVER_BEGIN -->
  :end-before: <!-- SPHINX_PROPPANT_COUPLED_SOLVER_END -->

References
=================================

- R. D. Barree & M. W. Conway. "Experimental and numerical modeling of convective proppant transport", JPT. Journal of petroleum technology, 47(3):216-222, 1995.

- M. A. Biot & W. L. Medlin. "Theory of Sand Transport in Thin Fluids", Paper presented at the SPE Annual Technical Conference and Exhibition, Las Vegas, NV, 1985.

- X. Hu, K. Wu, X. Song, W. Yu, J. Tang, G. Li, & Z. Shen. "A new model for simulating particle transport in a low-viscosity fluid for fluid-driven fracturing", AIChE J. 64 (9), 35423552, 2018.

- R. G. Keck, W. L. Nehmer, & G. S. Strumolo. "A new method for predicting friction pressures and rheology of proppant-laden fracturing fluids", SPE Prod. Eng., 7(1):21-28, 1992.

- M. McClure. "Bed load proppant transport during slickwater hydraulic fracturing: insights from comparisons between published laboratory data and correlations for sediment and pipeline slurry transport", J. Pet. Sci. Eng. 161 (2), 599610, 2018.

- M. C. Miller, I. N.  McCave, & P. D. Komar. "Threshold of sediment motion under unidirectional currents", Sedimentology 24 (4), 507527, 1977.

- P. L. Wiberg &  J. D. Smith. "Model for calculating bed load transport of sediment", J. Hydraul. Eng. 115 (1), 101123, 1989.
