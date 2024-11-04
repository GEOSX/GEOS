# ------------------------------------------------------------------------------------------------------------
# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 Total, S.A
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# Copyright (c) 2019-     INRIA project-team Makutu 
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
# ------------------------------------------------------------------------------------------------------------

import numpy as np
from mpi4py import MPI

from geos.pygeos_tools.utilities.input import XML
from geos.pygeos_tools.utilities.solvers import ReservoirSolver

from darts.models.darts_model import DartsModel
from darts.physics.super.physics import Compositional
from darts.physics.super.property_container import PropertyContainer
from darts.physics.properties.flash import ConstantK
from darts.physics.properties.basic import ConstFunc, PhaseRelPerm
from darts.physics.properties.density import DensityBasic

class Model(DartsModel):
    def __init__(self, n_points=50):
        # Call base class constructor
        super().__init__()
        self.n_obl_points = n_points
        self.set_physics()

    def set_physics(self):
        """Physical properties"""
        self.zero = 1e-8
        # Create property containers:
        components = ['CO2', 'C1', 'H2O']
        phases = ['gas', 'oil']
        thermal = 0
        Mw = [44.01, 16.04, 18.015]

        property_container = PropertyContainer(phases_name=phases, components_name=components,
                                               Mw=Mw, min_z=self.zero / 10, temperature=1.)

        """ properties correlations """
        property_container.flash_ev = ConstantK(len(components), [4, 2, 1e-1], self.zero)
        property_container.density_ev = dict([('gas', DensityBasic(compr=1e-3, dens0=200)),
                                              ('oil', DensityBasic(compr=1e-5, dens0=600))])
        property_container.viscosity_ev = dict([('gas', ConstFunc(0.05)),
                                                ('oil', ConstFunc(0.5))])
        property_container.rel_perm_ev = dict([('gas', PhaseRelPerm("gas")),
                                               ('oil', PhaseRelPerm("oil"))])

        """ Activate physics """
        self.physics = Compositional(components, phases, self.timer,
                                     n_points=self.n_obl_points, min_p=1, max_p=300, min_z=self.zero/10, max_z=1-self.zero/10)
        self.physics.add_property_region(property_container)
        self.engine = self.physics.init_physics(platform='cpu')
        return

def run_darts_model(xml_name: str, darts_model=None):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    xml = XML(xml_name)

    solver = ReservoirSolver()
    solver.initialize(rank=rank, xml=xml)
    solver.initialDt = solver.geosx.get_wrapper("Events/solver_1/forceDt").value()[0]
    solver.dt = solver.initialDt

    functions = solver.geosx.get_group("/Functions").groups()
    for func in functions:
        if hasattr(func, 'setAxes') and darts_model is not None:
          func.setAxes( darts_model.physics.n_vars, 
                        darts_model.physics.n_ops, 
                        list(darts_model.physics.axes_min), 
                        list(darts_model.physics.axes_max), 
                        list(darts_model.physics.n_axes_points) )
          func.setEvaluateFunction(darts_model.physics.reservoir_operators[0].evaluate)
          print("Adaptive OBL interpolator is configured.")

    if rank == 0:
       print()
       print("===============================================================")
       print(" Running pygeosx with solver:", solver.solver)
       print(" simulation duration: ", solver.maxTime)
       print(" initial timestep: ", solver.initialDt)
       print("===============================================================")
    solver.applyInitialConditions()

    time = 0
    cycle = 0
   
    solver.outputVtk(time)
    while time < solver.maxTime:
        if rank == 0:
          if solver.initialDt is not None:
              print(f"time = {time:.3f}s, dt = {solver.initialDt:.4f}, iter = {cycle+1}")
          else:
              print(f"time = {time:.3f}s, dt is None, iter = {cycle+1}")
        solver.execute(time)
        solver.updateTimeStep()
        solver.outputVtk(time)
        time += solver.initialDt
        cycle += 1
    solver.cleanup(time)
    
    comm.Barrier()


# run adaptive OBL
print("\n" + "="*30 + " RUNNING ADAPTIVE OBL " + "="*30 + "\n")
darts_model = Model()
run_darts_model(xml_name="input_file_adaptive.xml", darts_model=darts_model)

# run static OBL
print("\n" + "="*30 + " RUNNING STATIC OBL " + "="*30 + "\n")
run_darts_model(xml_name="input_file_static.xml")
