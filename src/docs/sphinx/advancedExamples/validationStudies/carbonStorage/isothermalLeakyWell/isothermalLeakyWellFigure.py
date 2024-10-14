import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py

#from matplotlib.ticker import AutoMinorLocator


def main():

    # File path
    hdf5FilePath = 'phaseOutfluxHistory.hdf5'

    # Read HDF5
    hf = h5py.File(hdf5FilePath, 'r')
    time = hf.get('phaseOutflux Time')
    phaseOutflux = hf.get('phaseOutflux')

    time = np.asarray(time)
    phaseOutflux = np.asarray(phaseOutflux)

    # Conversions
    inKilograms = 479.0
    inDays = 1.0 / 86400.0

    # Data from the benchmark
    timeStepSize = 10000
    injectionRate = 8.87

    eclipse_time = np.array([
        0, 1.389076824194035, 9.674909509150552, 33.72052704576976, 41.7991610567978, 44.5308007171611,
        48.681117012279685, 52.84834748486179, 61.195494063123704, 82.05278914786373, 114.02058455397315,
        166.8033896011637, 223.73862521565576, 337.562582456615, 477.71979973613884, 624.7906870538886,
        789.8772876425018, 914.7250600453298, 999.3424613511044
    ])
    eclipse_leakage = np.array([
        0, -0.00038109756097565173, 0.0060975609756097615, 0.171875, 0.21570121951219512, 0.22332317073170732,
        0.22522865853658536, 0.2240853658536585, 0.2195121951219512, 0.20998475609756095, 0.19778963414634143,
        0.18368902439024387, 0.17111280487804875, 0.15434451219512196, 0.1410060975609756, 0.13147865853658536,
        0.12461890243902438, 0.1208079268292683, 0.11852134146341461
    ])

    tough2_time = np.array([
        0.0021201488344502195, 11.162583613369634, 16.906066805889736, 22.73859624945672, 28.524482418665798,
        40.079293566408325, 45.68072678702043, 56.7987872748667, 72.06385888289353, 119.2668525330478,
        173.42817465786098, 294.31906119809605, 424.97959356746844, 632.1393361814, 808.737133346761, 915.8131300817319,
        999.2494673126054
    ])
    tough2_leakage = np.array([
        0, 0.006093307750204058, 0.03848441160569471, 0.08688157908685187, 0.12689461800218374, 0.20387192180891095,
        0.21072954321393358, 0.20920091590429646, 0.20309753744712877, 0.18783564606235353, 0.1733332979975194,
        0.15347015360478303, 0.13970084699945934, 0.12666458185364607, 0.12011862232728732, 0.11702903543828769,
        0.11471065268781863
    ])

    # Plot HDF5 content
    plt.plot(eclipse_time[:], eclipse_leakage[:], 'r--', label="ECLIPSE")
    plt.plot(tough2_time[:], tough2_leakage[:], 'b--', label="TOUGH2")
    plt.plot(time[:, 0] * inDays,
             100 * abs(phaseOutflux[:, 14, 0] / timeStepSize) * inKilograms / injectionRate,
             'k-',
             label="GEOSX")

    values = 100 * abs(phaseOutflux[:, 14, 0] / timeStepSize) * inKilograms / injectionRate
    argmax = np.argmax(values)

    plt.xlim(0, 1000)
    plt.ylim(0, 0.25)
    #ax = plt.axes()
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())

    plt.grid()
    plt.xlabel('Time [days]', fontweight='bold')
    plt.ylabel('Leakage value [%]', fontweight='bold')
    plt.legend(loc="upper right")
    plt.show()


if __name__ == "__main__":
    main()
