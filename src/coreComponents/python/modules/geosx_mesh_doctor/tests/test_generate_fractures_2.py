from dataclasses import dataclass

from typing import (
    Tuple,
    Iterable,
    Iterator,
    Dict,
    Mapping,
    FrozenSet,
    List,
    Set,
    Sequence,
    Collection,
)

import numpy

import pytest

from vtkmodules.vtkCommonDataModel import (
    vtkUnstructuredGrid,
)


from vtkmodules.util.numpy_support import (
    numpy_to_vtk,
    vtk_to_numpy,
)
import sys
sys.path.append("..")

# import sys
# sys.path.insert(0, "/Users/j0436735/CLionProjects/GEOS/src/coreComponents/python/modules/geosx_mesh_doctor")
# sys.path.insert(0, "/Users/j0436735/CLionProjects/GEOS/src/coreComponents/python/modules/geosx_mesh_doctor/tests")

from checks.generate_cube import build_rectilinear_blocks_mesh, XYZ
from checks.generate_fractures_2 import __split_mesh_on_fracture, Options
# from ..checks.generate_cube import build_rectilinear_blocks_mesh, XYZ
# from ..checks.generate_fractures_2 import __split_mesh_on_fracture, Options


@dataclass(frozen=True)
class TestCase:
    input_mesh: vtkUnstructuredGrid
    options: Options
    collocated_nodes: Tuple[Tuple[int, ...], ...]
    result: Tuple[int, ...]


def __build_test_case(xs: Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray],
                      attribute: Iterable[int],
                      field_values: Iterable[int] = None):
    xyz = XYZ(*xs)

    mesh: vtkUnstructuredGrid = build_rectilinear_blocks_mesh((xyz, ))

    ref = numpy.array(attribute, dtype=int)
    assert len(ref) == mesh.GetNumberOfCells()
    attr = numpy_to_vtk(ref)
    attr.SetName("attribute")
    mesh.GetCellData().AddArray(attr)

    if field_values is None:
        fv = frozenset(attribute)
    else:
        fv = field_values

    options = Options(policy="field",
                      field="attribute",
                      field_type="cells",
                      field_values=fv,
                      vtk_output=None,
                      vtk_fracture_output=None,
                      split_on_domain_boundary=True)
    return mesh, options


class Incrementor:
    def __init__(self, start):
        self.__val = start

    def next(self, num: int) -> Iterable[int]:
        self.__val += num
        return range(self.__val - num, self.__val)


def __generate_test_data() -> Iterator[Tuple[vtkUnstructuredGrid, Options]]:
    two_nodes = numpy.arange(2, dtype=float)
    three_nodes = numpy.arange(3, dtype=float)
    four_nodes = numpy.arange(4, dtype=float)

    # Split in 2
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), (0, 1, 0, 1, 0, 1, 0, 1))
    yield TestCase(input_mesh=mesh, options=options,
                   collocated_nodes=tuple(map(lambda i: (1 + 3 * i, 27 + i), range(9))),
                   result=(9 * 4, 8, 9, 4))

    # Split in 3
    inc = Incrementor(27)
    collocated_nodes = (
        (1, *inc.next(1)),
        (3, *inc.next(1)),
        (4, *inc.next(2)),
        (7, *inc.next(1)),
        (1 + 9, *inc.next(1)),
        (3 + 9, *inc.next(1)),
        (4 + 9, *inc.next(2)),
        (7 + 9, *inc.next(1)),
        (1 + 18, *inc.next(1)),
        (3 + 18, *inc.next(1)),
        (4 + 18, *inc.next(2)),
        (7 + 18, *inc.next(1)),
    )
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), (0, 1, 2, 1, 0, 1, 2, 1))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(9 * 4 + 6, 8, 12, 6))

    # Split in 8
    inc = Incrementor(27)
    collocated_nodes = (
        (1, *inc.next(1)),
        (3, *inc.next(1)),
        (4, *inc.next(3)),
        (5, *inc.next(1)),
        (7, *inc.next(1)),
        (0 + 9, *inc.next(1)),
        (1 + 9, *inc.next(3)),
        (2 + 9, *inc.next(1)),
        (3 + 9, *inc.next(3)),
        (4 + 9, *inc.next(7)),
        (5 + 9, *inc.next(3)),
        (6 + 9, *inc.next(1)),
        (7 + 9, *inc.next(3)),
        (8 + 9, *inc.next(1)),
        (1 + 18, *inc.next(1)),
        (3 + 18, *inc.next(1)),
        (4 + 18, *inc.next(3)),
        (5 + 18, *inc.next(1)),
        (7 + 18, *inc.next(1)),
    )
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), range(8))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(8 * 8, 8, 3 * 3 * 3 - 8, 12))

    # Straight notch
    inc = Incrementor(27)
    collocated_nodes = (
        (1, *inc.next(1)),
        (4,),
        (1 + 9, *inc.next(1)),
        (4 + 9,),
        (1 + 18, *inc.next(1)),
        (4 + 18,),
    )
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), (0, 1, 2, 2, 0, 1, 2, 2), field_values=(0, 1))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(3 * 3 * 3 + 3, 8, 6, 2))

    # L-shaped notch
    inc = Incrementor(27)
    collocated_nodes = (
        (1, *inc.next(1)),
        (4, *inc.next(1)),
        (7, *inc.next(1)),
        (1 + 9, *inc.next(1)),
        (4 + 9,),
        (7 + 9,),
        (1 + 18, *inc.next(1)),
        (4 + 18,),
    )
    mesh, options = __build_test_case((three_nodes, three_nodes, three_nodes), (0, 1, 0, 1, 0, 1, 2, 2), field_values=(0, 1))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(3 * 3 * 3 + 5, 8, 8, 3))

    # 3x1x1 split
    inc = Incrementor(2 * 2 * 4)
    collocated_nodes = (
        (1, *inc.next(1)),
        (2, *inc.next(1)),
        (5, *inc.next(1)),
        (6, *inc.next(1)),
        (1 + 8, *inc.next(1)),
        (2 + 8, *inc.next(1)),
        (5 + 8, *inc.next(1)),
        (6 + 8, *inc.next(1)),
    )
    mesh, options = __build_test_case((four_nodes, two_nodes, two_nodes), (0, 1, 2))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(6 * 4, 3, 2 * 4, 2))

    # No duplication one quad fracture
    collocated_nodes = (
        (4 + 12,),
        (7 + 12,),
        (4 + 12 * 2,),
        (7 + 12 * 2,),
    )
    mesh, options = __build_test_case((three_nodes, four_nodes, four_nodes), [0, ] * 8 + [1, 2] + [0, ] * 8, field_values=(1, 2))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(3 * 4 * 4, 2 * 3 * 3, 4, 1))

    # Fracture on a corner
    inc = Incrementor(3 * 4 * 4)
    collocated_nodes = (
        (1 + 12,),
        (4 + 12,),
        (7 + 12,),
        (1 + 12 * 2, *inc.next(1)),
        (4 + 12 * 2, *inc.next(1)),
        (7 + 12 * 2,),
        (1 + 12 * 3, *inc.next(1)),
        (4 + 12 * 3, *inc.next(1)),
        (7 + 12 * 3,),
    )
    mesh, options = __build_test_case((three_nodes, four_nodes, four_nodes), [0, ] * 6 + [1, 2, 1, 2, 0, 0, 1, 2, 1, 2, 0, 0], field_values=(1, 2))
    yield TestCase(input_mesh=mesh, options=options, collocated_nodes=collocated_nodes,
                   result=(3 * 4 * 4 + 4, 2 * 3 * 3, 9, 4))


def __format_collocated_nodes(fracture_mesh: vtkUnstructuredGrid):
    collocated_nodes: numpy.ndarray = vtk_to_numpy(fracture_mesh.GetPointData().GetArray("collocated_nodes"))
    if len(collocated_nodes.shape) == 1:
        collocated_nodes: numpy.ndarray = collocated_nodes.reshape((collocated_nodes.shape[0], 1))
    return tuple(sorted(map(lambda bucket: tuple(sorted(filter(lambda i: i != -1, bucket))), collocated_nodes)))  # TODO why the wrapping sorted?


@pytest.mark.parametrize("expected", __generate_test_data())
def test_generate_fracture(expected: TestCase):
    main_mesh, fracture_mesh = __split_mesh_on_fracture(expected.input_mesh, expected.options)
    assert main_mesh.GetNumberOfPoints() == expected.result[0]
    assert main_mesh.GetNumberOfCells() == expected.result[1]
    assert fracture_mesh.GetNumberOfPoints() == expected.result[2]
    assert fracture_mesh.GetNumberOfCells() == expected.result[3]

    res = __format_collocated_nodes(fracture_mesh)
    assert res == expected.collocated_nodes
    assert len(res) == expected.result[2]
