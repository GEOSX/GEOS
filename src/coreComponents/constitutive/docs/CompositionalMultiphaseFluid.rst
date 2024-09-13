.. _CompositionalMultiphaseFluid:

############################################
Compositional multiphase fluid model
############################################

Overview
=========================

This model represents a full composition description of a multiphase multicomponent fluid.
Phase behavior is modeled by an equation of state (EOS) and partitioning of components into
phases is computed based on instantaneous chemical equilibrium via a two-phase flash.
Each component (species) is characterized by molar weight and critical properties that
serve as input parameters for the EOS.
See `Petrowiki`_ for more information.

In this model the fluid is described by :math:`N_c` components with :math:`z_c` being the total
mole fraction of component :math:`c`. The fluid can partition into a liquid phase, denoted :math:`\ell`,
and a vapor phase denoted by :math:`v`. Therefore, by taking into account the molar phase component
fractions, (which is the fraction of the molar mass of phase :math:`p` represented by component 
:math:`c`), the following partition matrix establishes the component distribution within the two
phases:

.. math::
    \begin{bmatrix}
    x_{1} & x_{2} & x_{3} & \cdots & x_{N_c} \\
    y_{1} & y_{2} & y_{3} & \cdots & y_{N_c} \\
    \end{bmatrix}

where :math:`x_c` is the mole fraction of component :math:`c` in the liquid phase and :math:`y_c`
is the mole fraction of component :math:`c` in the vapor phase.

The fluid properties are updated through the following steps:

1) The phase fractions (:math:`\nu_p`) and phase component fractions (:math:`x_c` and :math:`y_c`) are
computed as a function of pressure (:math:`p`), temperature (:math:`T`) and total component fractions
(:math:`z_c`).

2) The phase densities (:math:`\rho_p`) and phase viscosities (:math:`\mu_p`) are computed as a function
of pressure, temperature and the updated phase component fractions.

After calculating the phase fractions, phase component fractions, phase densities, phase viscosities,
and their derivatives with respect to pressure, temperature, and component fractions, the
:ref:`CompositionalMultiphaseFlow` then moves on to assembling the accumulation and flux terms.

Step 1: Computation of the phase fractions and phase component fractions (flash)
================================================================================
Stability test
-------------------------------
The first step is to determine if the provided mixture with total molar fractions :math:`z_c` is stable
as a single phase at the current pressure :math:`p` and temperature :math:`T`. However, this can only
be confirmed through stability testing.

The stability of a mixture is traditionally assessed using the Tangent Plane Distance (TPD) criterion
developed by Michelsen (1982a). This criterion states that a phase with composition :math:`z` is stable
at a specified pressure :math:`p` and temperature :math:`T` if and only if 

.. math::
  g(y) = \sum_{i=1}^{N_c}y_i\left(\ln y_i + \ln \phi_i(y) - \ln z_i - \ln \phi_i(z) \right) \geq 0

for any permissible trial composition :math:`y`, where :math:`\phi_i` denotes the fugacity
coefficient of component :math:`i`. 

To determine stability of the mixture this testing in initiated from a basic starting point, based on
Wilson K-values, to get both a lighter and a heavier trial mixture. The two trial mixtures are
calculated as :math:`y_i = z_i/K_i` and :math:`y_i = z_iK_i` where :math:`K_i` are defined by

.. math::
  K_i = \frac{P_{ci}}{p}\exp\left( 5.37( 1 + \omega_i ) \left(1-\frac{T_{ci}}{T}\right)\right)
  
where :math:`P_{ci}` and :math:`T_{ci}` are respectively, the critical pressure and temperature of
component :math:`i` and :math:`\omega_i` is the accentric factor of component :math:`i`.

The stability problem is solved by observing that a necessary condition is that :math:`g(y)` must
be non-negative at all its stationary points. The stationarity criterion can be expressed as

.. math::
  \ln y_i + \ln \phi_i(y) - h_i = k \hspace{1cm} i=1,2,3,\ldots,N_c

where :math:`h_i = \ln z_i + \ln \phi_i(z)` is a constant parameter dependent on the feed composition
:math:`z` and :math:`k` is an undetermined constant. This constant can be further incorporated into
the equation by defining the unnormalized trial phase moles :math:`Y_i` as

.. math::
  Y_i = \exp(-k)y_i

which reduces the stationarity criterion to

.. math::
  \ln Y_i + \ln \phi_i(y) - h_i = 0

with the mole fractions :math:`y_i` related to the trial phase moles :math:`Y_i` by

.. math::
  y_i = Y_i/\sum_jY_j

With the two starting mixtures, the stationarity condition is solved using successive substitution to
determine the stationary points. If both initial states converge to a solution which has :math:`g(y)\geq 0`
then the mixture is deemed to be stable, otherwise it is deemed unstable.

Phase labeling
----------------------------------------
Once it is confirmed that the fluid with composition :math:`z` is stable as a single phase at the current
pressure and temperature, it must be labeled as either 'liquid' or 'vapor'. This is necessary only to apply
the correct relative permeability function for calculating the phase's flow properties. The properties of the
fluid (density, viscosity) are unchanged by the assignment of the label.

Determining the mixture's true critical point is the most rigorous method for labeling. It is however expensive
and may not always be necessary. As such, a simple correlation for pseudo-critical temperature is used and this
is expected to be sufficiently accurate for correct phase labeling, except under some specific conditions.

The Li-correlation is a weighted average of the component critical temperatures and is used to determine the label
applied to the mixture. The Li pseudo-critical temperature is calcaulated as

.. math::
  T_{cp} = \frac{\sum_{i=1}^{N_c}T_{ci}V_{ci}z_{i}}{\sum_{i=1}^{N_c}V_{ci}z_{i}}

where :math:`V_{ci}` and :math:`T_{ci}` are respectively the critical volume and temperature of component
:math:`i`. This is compared to the current temperature :math:`T` such that if :math:`T_{cp}<T` then the mixture
is labeled as vapor and as liquid otherwise.

Negative two-phase flash
------------------------
When a cell is identified as having an unstable mixture, it is necessary to determine the amounts in the liquid
and vapor phases through phase splitting. This phase split is calculated by ensuring that the two phases are
in thermodynamic equilibrium. For a system to be in thermodynamic equilibrium, the fugacities of each component
in both the liquid and vapor phases must be equal:

.. math::
    \phi_{iL} = \phi_{iV} \hspace{1cm} i=1,2,3,\ldots,N_c

where :math:`\phi_{iL}` is the fugacity of component :math:`i` in the liquid phase and :math:`\phi_{iL}` is the
fugacity of component :math:`i` in the vapor phase.

Fugacities are functions of temperature, pressure, and composition:

.. math::
    \phi_{iL} = \phi_{iL}(T, p, x_i) \hspace{1cm} i=1,2,3,\ldots,N_c

and 

.. math::
    \phi_{iV} = \phi_{iV}(T, p, y_i) \hspace{1cm} i=1,2,3,\ldots,N_c

and are calculated directly from an equation of state.

Equilibrium constants, also known as K-values, are defined for each component as:

.. math::
    K_i = \frac{y_i}{x_i}

where :math:`x_i` is the mole fraction of component :math:`i` in the liquid phase and :math:`y_i` is the
mole fraction of component :math:`i` in the vapor phase. If we denote :math:`V` as the mole fraction
of the vapor phase, the material balance indicates that the mole fractions of each component in the liquid
and vapor phases are given by:

.. math::
    x_i = \frac{z_i}{1 + (K_i - 1)V}

and

.. math::
    y_i = \frac{K_i z_i}{1 + (K_i - 1)V}

The value of :math:`V` corresponding to a given set of K-values is determined by solving the
so called Rachford and-Rice equation:

.. math::
    F(V) = \sum_{i=1}^{N_c} \left(x_i - y_i\right) = \sum_{i=1}^{N_c} \frac{z_i(1 - K_i)}{1 + (K_i - 1)V} = 0

The flash calculation process is as follows:

#. Once the mixture is confirmed to be stable, an initial set of K-values is chosen, typically using Wilson's formula.

#. Given :math:`z_i` and :math:`K_i`, the Rachford-Rice equation is solved to determine the molar fraction of vapor,  :math:`V`. This is initially solved using successive substitution, followed by Newton iterations once the residual is sufficiently reduced.

#. After  :math:`V` is calculated, the corresponding liquid and vapor mole fractions, :math:`x_i` and :math:`y_i`, are computed.

#. These phase compositions are then used to calculate the component fugacities :math:`\phi_{iL}` and :math:`\phi_{iV}` in the liquid and vapor phases using the equation of state.

#. Convergence is reached when the fugacities are equal for all components. The convergence criterion is defined as:

   .. math::
       \sum_{i=1}^{N_c} \left( \phi_{iL} - \phi_{iV} \right)^2 < \varepsilon
   
   where :math:`\varepsilon` is the convergence tolerance.

#. If convergence is not achieved, successive substitution is used to update the set of K-values for the next iteration. The new K-values at iteration  :math:`t+1` are given by:

   .. math::
       K_i^{(t+1)} = K_i^{(t)} \frac{\phi_{iL}}{\phi_{iV}}

Parameters
=========================

The model represented by ``<CompositionalMultiphaseFluid>`` node in the input.
Under the hood this is a wrapper around ``PVTPackage`` library, which is included as a submodule.
In order to use the model, GEOS must be built with ``-DENABLE_PVTPACKAGE=ON`` (default).

The following attributes are supported:

.. include:: ../../../coreComponents/schema/docs/CompositionalMultiphaseFluid.rst

Supported phase names are:

===== ===========
Value Comment
===== ===========
oil   Oil phase
gas   Gas phase
water Water phase
===== ===========

Supported Equation of State types:

===== =======================
Value Comment
===== =======================
PR    Peng-Robinson EOS
SRK   Soave-Redlich-Kwong EOS
===== =======================

Example
=========================

.. code-block:: xml

  <Constitutive>
    <CompositionalMultiphaseFluid name="fluid1"
                                  phaseNames="{ oil, gas }"
                                  equationsOfState="{ PR, PR }"
                                  componentNames="{ N2, C10, C20, H2O }"
                                  componentCriticalPressure="{ 34e5, 25.3e5, 14.6e5, 220.5e5 }"
                                  componentCriticalTemperature="{ 126.2, 622.0, 782.0, 647.0 }"
                                  componentAcentricFactor="{ 0.04, 0.443, 0.816, 0.344 }"
                                  componentMolarWeight="{ 28e-3, 134e-3, 275e-3, 18e-3 }"
                                  componentVolumeShift="{ 0, 0, 0, 0 }"
                                  componentBinaryCoeff="{ { 0, 0, 0, 0 },
                                                        { 0, 0, 0, 0 },
                                                        { 0, 0, 0, 0 },
                                                        { 0, 0, 0, 0 } }"/>
  </Constitutive>

References
==========

- M. L. Michelsen, `The Isothermal Flash Problem. Part I. Stability.
  <https://doi.org/10.1016/0378-3812(82)85001-2>`__, Fluid Phase Equilibria,
  vol. 9.1, pp. 1-19, 1982a.

- M. L. Michelsen, `The Isothermal Flash Problem. Part II. Phase-Split Calculation.
  <https://doi.org/10.1016/0378-3812(82)85002-4>`__, Fluid Phase Equilibria,
  vol. 9.1, pp. 21-40, 1982b.

.. _Petrowiki: https://petrowiki.spe.org/Phase_behavior_in_reservoir_simulation#Equation-of-state_models
