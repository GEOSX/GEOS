.. _DelftEgg:

=========================
Model: Delft Egg
=========================
The `Delft-Egg  <https://link.springer.com/chapter/10.1007%2F978-94-011-1046-4_10>`__  plasticity model uses a generalization of the Modified Cam-Clay yield surface, defined as 

.. math::
  f = q^2 - M^2 \left[ \alpha^2 \, p \left(\frac{2 \alpha}{\alpha+1} p_c -p \right) - \frac{\alpha^2 (\alpha-1)}{\alpha+1} p_c^2 \right] = 0 \quad ( \text{for }p_c > \frac{\alpha}{\alpha+1} )
  
.. math::
  f = q^2 - M^2 \, p \left(\frac{2 \alpha}{\alpha+1} p_c -p \right) = 0 \quad (\text{for } p_c \leq \frac{\alpha}{\alpha+1} )

where :math:`\alpha \geq 1` is the shape parameter. For :math:`\alpha = 1`, this model leads to a Modified Cam-Clay (MCC) type model with an ellipsoidal yield surface. For :math:`\alpha > 1`, an egg-shaped yield surface is obtained.  The additional parameter makes it easier to fit the cap behavior of a broader range of soils and rocks. 

Because Delft-Egg is frequently used for hard rocks, GEOS uses a linear model for the elastic response, rather than the hyper-elastic model used for MCC.  This is a slight deviation from the original formulation proposed in the reference above.  For reservoir applications, the ability to use a simpler linear model was a frequent user request.

Parameters
~~~~~~~~~~~~~~~~~~~~

The supported attributes will be documented soon.


Example
~~~~~~~~~~~~~~~~~~~

.. code-block:: xml

  <Constitutive>
     <DelftEgg
      name="DE"
      defaultDensity="2700"
      defaultBulkModulus="10.0e9"
      defaultShearModulus="6.0e9"     
      defaultPreConsolidationPressure="-20.0e6"
      defaultShapeParameter="6.5"
      defaultCslSlope="1.2"
      defaultVirginCompressionIndex="0.005"
      defaultRecompressionIndex="0.001"/>
  </Constitutive>
