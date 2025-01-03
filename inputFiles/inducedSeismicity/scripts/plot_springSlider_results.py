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
    var = data[f'{var_name}']
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
    
    filePath = f"{normalized_dir}/springSlider.hdf5"
    if not os.path.exists(filePath):
        print(f"Error: {filePath} not found.")
        exit(1)    
    
    time, slipRate = getDataFromHDF5( filePath, "slipRate" )
    time, stateVariable = getDataFromHDF5( filePath, "stateVariable" )

    # Convert time to years
    time_in_years = time / (365 * 24 * 3600)  # Assuming time is in seconds

    # Plotting
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

    # Plot pressure on the left y-axis
    ax1.set_xlabel('Time (years)')
    ax1.set_ylabel('slip rate (Pa)', color='tab:blue')
    ax1.plot(time_in_years, slipRate, label="slip rate", color='tab:blue')
    ax1.set_yscale('log')
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    # Set x-axis limits to 0 to 2 years
    ax1.set_xlim(0, np.max(time_in_years))
  
    # Plot stateVariable on the right y-axis (log scale)
    ax2.set_ylabel('state variable', color='tab:red')
    ax2.plot(time_in_years, stateVariable, label="state variable", color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')
    ax2.set_xlim(0, np.max(time_in_years))

    # Add grid and title
    fig.suptitle('Spring slider solution', fontsize=16)
    plt.tight_layout()
    plt.grid()

    # Show plot
    plt.show()