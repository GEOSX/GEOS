.. _dfits:


####################################################
Diagnostic Fracture-Injection Tests
####################################################


------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

 
The parameters used in the simulation are summarized in the following table.

+------------------+-------------------------+--------------------+--------------------+
| Symbol           | Parameter               | Unit               | Value              |
+==================+=========================+====================+====================+
| :math:`K`        | Bulk Modulus            | [GPa]              | 4.11               |
+------------------+-------------------------+--------------------+--------------------+
| :math:`G`        | Shear Modulus           | [GPa]              | 1.2                |
+------------------+-------------------------+--------------------+--------------------+
| :math:`K_{Ic}`   | Rock Toughness          | [MPa.m\ :sup:`1/2`]| 1.2                |
+------------------+-------------------------+--------------------+--------------------+
| :math:`\mu`      | Fluid Viscosity         | [Pa.s]             | 97.7               |
+------------------+-------------------------+--------------------+--------------------+
| :math:`Q_0`      | Injection Rate          | [m\ :sup:`3`/s]    | 73.2x10\ :sup:`-9` |
+------------------+-------------------------+--------------------+--------------------+
| :math:`t_{inj}`  | Injection Time          | [s]                | 100                |
+------------------+-------------------------+--------------------+--------------------+
| :math:`h_f`      | Fracture Height         | [mm]               | 55                 |
+------------------+-------------------------+--------------------+--------------------+

**Input file**


.. code-block:: console

  ifixThis.xml

.. code-block:: console

  fixThis.xml

------------------------------------------------------------------
Results
------------------------------------------------------------------

Python scripts for post-processing and visualizing the simulation results are also prepared:

.. code-block:: console

  fixThis.py

.. code-block:: console

  fixThis.py


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
