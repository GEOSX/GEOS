from mpi4py import MPI
import pygeosx
from geosx_xml_tools.main import preprocess_parallel


def run_problem():
    """
    Run the GEOSX problem
    """
    # PYGEOSX_INITIALIZATION
    # Get the MPI rank
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Preprocess the xml file
    args = preprocess_parallel()

    # Initialize the code and set initial conditions
    problem = pygeosx.initialize(rank, args)
    pygeosx.apply_initial_conditions()
    while pygeosx.run() != pygeosx.COMPLETED:
        pass


if __name__ == '__main__':
    run_problem()
