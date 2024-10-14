import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py


def main():

    numWells = 4
    iplt = -1
    cmap = plt.get_cmap("tab10")

    # Loop over the four producers
    for iw in range(1, numWells + 1):
        # File path
        hdf5FilePath = 'wellRateHistory' + str(iw) + '.hdf5'

        # Read HDF5
        hf = h5py.File(hdf5FilePath, 'r')
        time = hf.get('Time')
        time = np.asarray(time)
        massRate = hf.get('wellElementMixtureConnectionRate')
        massRate = np.asarray(massRate)

        # Some comments about the computation of the volumetric rate here:
        # A proper oil rate constraint for individual wells is currently being implemented
        # In the meantime, the volume rate is (wrongly) computed by dividing
        # the total mass rate by the surface oil density

        # Conversions
        inCubicMeters = 1 / 848.9
        inDays = 1.0 / (24 * 3600)

        # Plot HDF5 content (here, the rate at the well)
        iplt += 1
        plt.plot(time[:, 0] * inDays,
                 abs(massRate[:, 0]) * inCubicMeters / inDays,
                 '-o',
                 color=cmap(iplt),
                 label='Producer #' + str(iw))

    plt.xlim(-1, 175)
    plt.ylim(0, 3800)
    plt.grid()
    plt.xlabel('time [days]')
    plt.ylabel('total rate [cubic meters per day]')
    plt.legend(bbox_to_anchor=(0.025, 0.975), loc='upper left', borderaxespad=0.)
    plt.show()


if __name__ == "__main__":
    main()
