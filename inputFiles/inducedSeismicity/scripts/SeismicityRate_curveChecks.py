
import argparse
import numpy as np
import xml.etree.ElementTree as ElementTree
import matplotlib.pyplot as plt

class analyticalSolution():
    def __init__( self, parameters ):
        self.cTau = 1e-3
        self.aSigma = 1e6
        self.t_a = 1e6/parameters['backgroundStressingRate']
        self.initial_normal_traction     = parameters['initialFaultNormalTraction']
        self.initial_fault_shear_traction = parameters['initialFaultShearTraction']

    def compute_shear_stress( self, time_n, dt ):
        return self.aSigma * np.log( self.cTau*( time_n + dt ) + 1 ) + self.initial_fault_shear_traction

    def compute_seismic_rate(self, time_n, dt, shear_stress):
        K = np.exp( shear_stress / self.aSigma - self.initial_fault_shear_traction / self.aSigma )
        denom = ( self.cTau*(time_n + dt - self.t_a ) + 1) * np.exp( (time_n+dt)  / self.t_a ) + self.cTau*self.t_a
        return K/denom

def getParametersFromXML( xmlFilePath ):
    tree = ElementTree.parse( xmlFilePath )
    root = tree.getroot()
    
    parameters = dict.fromkeys(["initialFaultNormalTraction", "initialFaultShearTraction", "backgroundStressingRate"])

    solver = tree.find('Solvers/SeismicityRate')
    parameters["backgroundStressingRate"]  = float(solver.get("backgroundStressingRate"))

    for field_spec in root.findall(".//FieldSpecification"):
        if field_spec.get('name') == 'initialShearTraction':
            parameters["initialFaultShearTraction"] = float(field_spec.get("scale")) 
        elif field_spec.get('name') == 'initialNormalTraction':
            parameters["initialFaultNormalTraction"]  = float(field_spec.get("scale")) 

    return parameters

def curve_check_solution(**kwargs):
     
    xmlFilePath = "./SeismicityRate_analytical_verification_smoke.xml"
    parameters = getParametersFromXML( xmlFilePath )

    analytical_solution = analyticalSolution( parameters )

    times = np.squeeze(kwargs['seismicityRate Time'][:])
    analytical_rates = []
    dt = 3600.0
    for time in times:
        tau = analytical_solution.compute_shear_stress(time, dt)
        analytical_rate = analytical_solution.compute_seismic_rate(time, dt, tau)
        analytical_rates.append(analytical_rate)

    return np.array(analytical_rates)   


def debug( xmlFilePath ):
    #-------- Extract info from XML
    parameters = getParametersFromXML( xmlFilePath )

    analytical_solution = analyticalSolution( parameters )

    time = 0.0
    dt = 3600.0
    final_time = 360000.0
    times = []
    analytical_rates = []
    tau_plot = []
    while time < final_time:
        tau = analytical_solution.compute_shear_stress(time, dt)
        analytical_rate = analytical_solution.compute_seismic_rate(time, dt, tau)
        time += dt
        times.append(time)
        analytical_rates.append(analytical_rate)
        tau_plot.append(tau)

    import csv

    # Write times to a CSV file without headers
    with open('shearStress_time.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        for time in times:
            writer.writerow([time])

    # Write tau_plot to a CSV file without headers
    with open('shearStress_values.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        for tau in tau_plot:
            writer.writerow([tau])

    import h5py            
    file_path = '/usr/workspace/cusini1/geosx/geosx_dev/GEOS_2/build-quartz-gcc-12-release/Output/seismicityRate.hdf5'
    with h5py.File(file_path, 'r') as file:
        # List all groups
        print("Keys: %s" % file.keys())

        # Get the data
        time = np.squeeze(file['seismicityRate Time'][:])
        seismicityRate = np.squeeze(file['seismicityRate'][:])
        print(time)
        print(seismicityRate)
    
    # Plot analytical (continuous line) and numerical (markers) aperture solution 
    fig, ax = plt.subplots(figsize=(16, 12), nrows=2, ncols=1)
    
    ax[0].plot(times, tau_plot)
    ax[0].set_xlabel('time [s]', weight="bold")
    ax[0].set_ylabel('shear stress', weight="bold")
    
    ax[1].plot(times, analytical_rates)
    ax[1].set_xlabel('time [s]', weight="bold")
    ax[1].set_ylabel('seismic rate', weight="bold")
    plt.savefig("seismicRate.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--xml-file', type=str, help='Path to XML file')
    args = parser.parse_args()
    debug( args.xml_file )
