import pygeosx
import sys
from mpi4py import MPI

# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def print_and_flush(msg, *args, **kwargs):
    print("Rank {}: {}".format(rank, msg), *args, **kwargs)
    sys.stdout.flush()


def print_with_indent(msg, indent):
    indent_str = " " * indent
    print(indent_str + msg.replace("\n", "\n" + indent_str))


def print_group(group, indent=0):
    print("{}{}".format(" " * indent, group))

    indent += 4
    print("{}wrappers:".format(" " * indent))

    for wrapper in group.wrappers():
        print("{}{}".format(" " * (indent + 4), wrapper))
        print_with_indent(str(wrapper.value(False)), indent + 8)

    print("{}groups:".format(" " * indent))

    for subgroup in group.groups():
        print_group(subgroup, indent + 4)


def callback(matrix, array):
    pass


def main():
    print_and_flush("In python before initialization.")

    problem = pygeosx.initialize(rank, sys.argv)

    try:
        problem.get_group("Solvers/lagsolve", None).register(callback)
    except AttributeError:
        pass

    curr_time = problem.get_wrapper("Events/time").value(False)
    print_and_flush("In python after initialization: current time = {}".format(curr_time[0]))

    pygeosx.apply_initial_conditions()
    curr_time = problem.get_wrapper("Events/time").value(False)
    print_and_flush("In python after applyInitialConditions: current time = {}".format(curr_time[0]))

    while pygeosx.run() != pygeosx.COMPLETED:
        curr_time = problem.get_wrapper("Events/time").value(True)
        print_and_flush("In python: current time = {}".format(curr_time[0]))
        curr_time[0] += 1e-6

    curr_time = problem.get_wrapper("Events/time").value(False)
    print_and_flush("In python after after the simulation has ended: current time = {}".format(curr_time[0]))

    # node_manager = problem.get_group("domain/MeshBodies/mesh1/Level0/nodeManager")
    # print_group(node_manager)
    # x = node_manager.get_wrapper("totalDisplacement")


if __name__ == "__main__":
    main()
