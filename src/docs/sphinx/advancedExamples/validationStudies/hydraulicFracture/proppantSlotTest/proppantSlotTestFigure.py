import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()
    outputDir = args.outputDir

    eTime, eLength, eArea = np.loadtxt(outputDir + "/experiment-data.txt", skiprows=0, unpack=True)
    esTime, esArea = np.loadtxt(outputDir + "/experiment-data2.txt", skiprows=0, unpack=True)

    dTime, dLength, dArea = np.loadtxt(outputDir + "/model-results.txt", skiprows=1, unpack=True)
    sTime, sArea = np.loadtxt(outputDir + "/model-results2.txt", skiprows=1, unpack=True)

    # offSize is the size of the bounday cell which should be excluded from the model domain
    offSize = 0.0127

    # slot length 4 feet = 1.2192 m
    totalLength = 1.2192

    #total cell number of the domain (excluding the bounday cells)
    totalCellNum = 96 * 24
    font = {'size': 12}
    matplotlib.rc('font', **font)

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))

    ax[0, 0].plot(eTime, eLength, 'ro', label='Experiment-30/50')
    ax[0, 0].plot(dTime, (dLength - offSize) / totalLength,
                  lw=5,
                  alpha=0.5,
                  color='lime',
                  label='GEOSX-30/50 mesh',
                  linestyle='-')
    ax[0, 0].set_ylabel('Normalized Bank Length', multialignment='center', weight="bold")
    ax[0, 0].set_xlabel('Time (s)', weight="bold")
    ax[0, 0].set_xlim(0, 40)
    ax[0, 0].set_ylim(0.0, 1.2)
    ax[0, 0].legend(frameon=False, loc='lower right', fontsize=12)
    ax[0, 0].text(.01, .92, '(a)', horizontalalignment='left', transform=ax[0, 0].transAxes, size=15, weight='bold')

    ax[0, 1].plot(eTime, eArea, 'ro', label='Experiment-30/50')
    ax[0, 1].plot(dTime, dArea / totalCellNum, lw=5, alpha=0.5, color='lime', label='GEOSX-30/50 mesh', linestyle='-')
    ax[0, 1].set_ylabel('Normalized Bank Area', multialignment='center', weight="bold")
    ax[0, 1].set_xlabel('Time (s)', weight="bold")
    ax[0, 1].set_xlim(0, 40)
    ax[0, 1].set_ylim(0.0, 0.5)
    ax[0, 1].legend(frameon=False, loc='lower right', fontsize=12)
    ax[0, 1].text(.01, .92, '(b)', horizontalalignment='left', transform=ax[0, 1].transAxes, size=15, weight='bold')

    ax[1, 0].plot(esTime, esArea, 'ro', label='Experiment-30/50')
    ax[1, 0].plot(sTime, sArea / totalCellNum, lw=5, alpha=0.5, color='lime', label='GEOSX-30/50 mesh', linestyle='-')
    ax[1, 0].set_ylabel('Normalized Suspended Proppant Area', multialignment='center', weight="bold")
    ax[1, 0].set_xlabel('Time (s)', weight="bold")
    ax[1, 0].set_xlim(0, 40)
    ax[1, 0].set_ylim(0.0, 0.6)
    ax[1, 0].legend(frameon=False, loc='upper right', fontsize=12)
    ax[1, 0].text(.01, .92, '(c)', horizontalalignment='left', transform=ax[1, 0].transAxes, size=15, weight='bold')

    ax[1, 1].plot(eTime, eArea + esArea, 'ro', label='Experiment-30/50')
    ax[1, 1].plot(dTime, (dArea + sArea) / totalCellNum,
                  lw=5,
                  alpha=0.5,
                  color='lime',
                  label='GEOSX-30/50 mesh',
                  linestyle='-')
    ax[1, 1].set_ylabel('Normalized Total Propped Area', multialignment='center', weight="bold")
    ax[1, 1].set_xlabel('Time (s)', weight="bold")
    ax[1, 1].set_xlim(0, 40)
    ax[1, 1].set_ylim(0.0, 1.0)
    ax[1, 1].legend(frameon=False, loc='upper right', fontsize=12)
    ax[1, 1].text(.01, .92, '(d)', horizontalalignment='left', transform=ax[1, 1].transAxes, size=15, weight='bold')

    plt.show()


if __name__ == "__main__":
    main()
