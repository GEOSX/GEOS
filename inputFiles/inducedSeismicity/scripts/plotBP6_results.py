import argparse
import numpy as np
from geos.hdf5_wrapper import hdf5_wrapper
import matplotlib.pyplot as plt
import os

def remove_padding(data):
    """Removes trailing zeros from a NumPy array."""
    nonzero_indices = np.nonzero(data)[0]
    if nonzero_indices.size == 0:  # If all elements are zero
        return np.array([]), np.array([])
    last_nonzero = nonzero_indices[-1]
    return data[:last_nonzero + 1]

def getDataFromHDF5( hdf5FilePath, var_name ):
    # Read HDF5
    data = hdf5_wrapper(f'{hdf5FilePath}').get_copy()
    var = data[f'{var_name} source']
    var = np.asarray(var)
    time = data[f'{var_name} Time']
    
    # Remove padding
    var = remove_padding(var)
    time = time[:len(var)]  # Match the length of the time array to the variable
    return time, var    


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
    
    time, pressure = getDataFromHDF5( filePath, "pressure" )
    time, slipRate = getDataFromHDF5( filePath, "slipRate" )

    # Convert time to years
    time_in_years = time / (365 * 24 * 3600)  # Assuming time is in seconds

    # Plotting
    fig, ax1 = plt.subplots()
    
    print(len(time_in_years))
    print(len(pressure)) 

    # Plot pressure on the left y-axis
    ax1.set_xlabel('Time (years)')
    ax1.set_ylabel('Pressure (Pa)', color='tab:blue')
    ax1.plot(time_in_years, pressure, label="Pressure", color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
  

    # Plot slipRate on the right y-axis (log scale)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Slip Rate (m/s)', color='tab:red')
    ax2.plot(time_in_years, slipRate, label="Slip Rate", color='tab:red')
    ax2.set_yscale('log')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    # Set x-axis limits to 0 to 2 years
    ax1.set_xlim(0, np.max(time_in_years))


    # Add grid and title
    plt.title("Pressure and Slip Rate vs Time")
    plt.grid()

    # Show plot
    plt.show()