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


    # Experiments (Rubin, 1983)
    H = 0.055
    Tinj = 3.7
    tlist = [3.7, 11.9, 12.55, 15.5, 16.05, 19.25, 24.30, 30.35, 38.50, 46.15, 53.85, 61.50, 68.90]
    Pblist = [0, 5.43, 5.26, 4.20, 4.01, 3.23, 2.51, 2.03, 1.61, 1.39, 1.22, 1.07, 0.97]
    lenlist = [24, 30, 31, 40, 43, 50, 60, 70, 81, 90, 97, 105, 112]
    widlist = [25, 69, 78, 133, 142, 185, 236, 294, 364, 420, 470, 514, 552]
    P57list = [0.371, 0.362, 0.365, 0.362, 0.359, 0.412, 0.886, 1.396, 1.279, 1.153, 1.034, 0.942, 0.861]
    P58list = [1.677, 4.145, 4.134, 3.507, 3.381, 2.827, 2.292, 1.894, 1.557, 1.345, 1.172, 1.050, 0.949]

    # Load GEOSX results
    GTime, GWellP, G58P, G57P, GAper, GArea = np.loadtxt(outputDir + "/model_results.txt", skiprows=1, unpack=True)
    GLength = GArea / H * 1000
    GTime = GTime + Tinj

    # Visulization
    N1 = 1
    fsize = 20
    msize = 15
    lw = 8
    lablelist = ['Experiment (Rubin, 1983)', 'GEOSX']
    fig, ax = plt.subplots(3, 2, figsize=(32, 18))
    cmap = plt.get_cmap("tab10")

    ax[0, 0].plot(tlist,
                  Pblist,
                  'o',
                  alpha=1.0,
                  color=cmap(0),
                  mec=cmap(0),
                  fillstyle='none',
                  markersize=msize,
                  mew=5,
                  label=lablelist[0])
    ax[0, 0].plot(GTime, GWellP / 1000000, lw=lw, alpha=0.5, color=cmap(2), label=lablelist[1])
    ax[0, 0].set_xlim([0, max(GTime)])
    ax[0, 0].set_ylim(0, 10)
    ax[0, 0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[0, 0].set_ylabel(r'Pressure @ Well (MPa)', size=fsize, weight="bold")
    ax[0, 0].legend(loc='upper right', fontsize=fsize)
    ax[0, 0].grid(True)
    ax[0, 0].xaxis.set_tick_params(labelsize=fsize)
    ax[0, 0].yaxis.set_tick_params(labelsize=fsize)

    ax[1, 0].plot(tlist,
                  P58list,
                  'o',
                  alpha=1.0,
                  color=cmap(0),
                  mec=cmap(0),
                  fillstyle='none',
                  markersize=msize,
                  mew=5,
                  label=lablelist[0])
    ax[1, 0].plot(GTime[13:], G58P[13:] / 1000000, lw=lw, alpha=0.5, color=cmap(2), label=lablelist[1])
    ax[1, 0].set_xlim([0, max(GTime)])
    ax[1, 0].set_ylim(0, 5)
    ax[1, 0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[1, 0].set_ylabel(r'Pressure @ Gage 58 (MPa)', size=fsize, weight="bold")
    ax[1, 0].grid(True)
    ax[1, 0].xaxis.set_tick_params(labelsize=fsize)
    ax[1, 0].yaxis.set_tick_params(labelsize=fsize)

    ax[2, 0].plot(tlist,
                  P57list,
                  'o',
                  alpha=1.0,
                  color=cmap(0),
                  mec=cmap(0),
                  fillstyle='none',
                  markersize=msize,
                  mew=5,
                  label=lablelist[0])
    ax[2, 0].plot(GTime, G57P / 1000000, lw=lw, alpha=0.5, color=cmap(2), label=lablelist[1])
    ax[2, 0].set_xlim([0, max(GTime)])
    ax[2, 0].set_ylim(0, 5)
    ax[2, 0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[2, 0].set_ylabel(r'Pressure @ Gage 57 (MPa)', size=fsize, weight="bold")
    ax[2, 0].grid(True)
    ax[2, 0].xaxis.set_tick_params(labelsize=fsize)
    ax[2, 0].yaxis.set_tick_params(labelsize=fsize)

    ax[0, 1].axis('off')

    ax[1, 1].plot(tlist,
                  lenlist,
                  'o',
                  alpha=1.0,
                  color=cmap(0),
                  mec=cmap(0),
                  fillstyle='none',
                  markersize=msize,
                  mew=5,
                  label=lablelist[0])
    ax[1, 1].plot(GTime, GLength, lw=lw, alpha=0.5, color=cmap(2), label=lablelist[1])
    ax[1, 1].set_xlim([0, max(GTime)])
    ax[1, 1].set_ylim(0, 150)
    ax[1, 1].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[1, 1].set_ylabel(r'Fracture Length (mm)', size=fsize, weight="bold")
    ax[1, 1].grid(True)
    ax[1, 1].xaxis.set_tick_params(labelsize=fsize)
    ax[1, 1].yaxis.set_tick_params(labelsize=fsize)

    ax[2, 1].plot(tlist,
                  widlist,
                  'o',
                  alpha=1.0,
                  color=cmap(0),
                  mec=cmap(0),
                  fillstyle='none',
                  markersize=msize,
                  mew=5,
                  label=lablelist[0])
    ax[2, 1].plot(GTime, GAper * 1e6, lw=lw, alpha=0.5, color=cmap(2), label=lablelist[1])
    ax[2, 1].set_xlim([0, max(GTime)])
    ax[2, 1].set_ylim(0, 800)
    ax[2, 1].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[2, 1].set_ylabel(r'Aperture @ LVDT (um)', size=fsize, weight="bold")
    ax[2, 1].grid(True)
    ax[2, 1].xaxis.set_tick_params(labelsize=fsize)
    ax[2, 1].yaxis.set_tick_params(labelsize=fsize)

    plt.show()


if __name__ == "__main__":
    main()
