import argparse
import numpy as np
from geos.hdf5_wrapper import hdf5_wrapper
import matplotlib.pyplot as plt
import os
from scipy.special import erfc

def remove_padding(data):
    """Removes trailing zeros from a NumPy array."""
    nonzero_indices = np.nonzero(data)[0]
    if nonzero_indices.size == 0:  # If all elements are zero
        return np.array([]), np.array([])
    last_nonzero = nonzero_indices[-1]
    return data[:last_nonzero + 1]

def getDataFromHDF5( hdf5FilePath, var_name, set_name):
    # Read HDF5
    data = hdf5_wrapper(f'{hdf5FilePath}').get_copy()
    var = data[f'{var_name} {set_name}']
    var = np.asarray(var)
    time = data[f'{var_name} Time']
    
    # Remove padding
    var = remove_padding(var)
    time = time[:len(var)]  # Match the length of the time array to the variable
    return time, var 

def analytical_pressure( time, x ):
    alpha = 0.1 # m^2 / s
    phi = 0.1
    beta = 1e-8 # Pa^-1
    t_off = 100 * 24 * 3600
    q0 = 1.25 * 1e-6 #m/s

    pressure = q0 / (beta * phi * np.sqrt(alpha)) * (G(time, x, alpha) - np.where(time > t_off, G(time - t_off, x, alpha), 0))
    return pressure 

def G( t, x, alpha ):
    A = np.abs(x) / np.sqrt( 4*alpha*t )
    B = -x**2 / (4 * alpha * t)
    return np.sqrt(t) * ( np.exp(B) / np.sqrt(np.pi)  - A * erfc( A ) )
     
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', type=str, help='Path to hdf5 file')
    args = parser.parse_args()

    normalized_dir = os.path.abspath(args.dir) 
    if not os.path.isdir(normalized_dir):
        print(f"Error: {normalized_dir} is not a valid directory.")
        exit(1)
    
    filePath = f"{normalized_dir}/BP6_DQ_S.hdf5"
    if not os.path.exists(filePath):
        print(f"Error: {filePath} not found.")
        exit(1)    
    
    # Plotting
    _, ax1 = plt.subplots()
    # _, ax2 = plt.subplots()
    
    positions_along_fault = [0., 500., 1500., 2500., 3500., 5000., 7500., -1500.]
    set_names = ["source", "receiver1", "receiver2", "receiver3", "receiver4", "receiver5", "receiver6", "receiver7"]
    # Plot pressure on the left y-axis
    for position, set_name in zip(positions_along_fault, set_names):
        time, pressure = getDataFromHDF5( filePath, "pressure" , set_name)
        time_in_years = time / (365 * 24 * 3600)  # Convert time to years, assuming time is in seconds
        pressure_analytical = analytical_pressure( time, position )
        ax1.plot(time_in_years, pressure, label=f"Pressure z = {position} m")
        ax1.plot(time_in_years, pressure_analytical, label=f"Pressure Analytical z = {position} m", linestyle='--')
        ax1.set_xlabel('Time (years)')
        ax1.set_ylabel('Pressure (Pa)', color='tab:blue')
        ax1.tick_params(axis='y', labelcolor='tab:blue')
        ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        

    # Set x-axis limits to 0 to 2 years
    ax1.set_xlim(0, np.max(time_in_years))
    # ax2.set_xlim(0, np.max(time_in_years))


    # Add grid and title
    plt.title("Pressure and Slip Rate vs Time")
    plt.grid()
    plt.tight_layout()

    # Show plot
    plt.show()