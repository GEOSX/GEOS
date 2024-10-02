import sys
import os
import argparse


sys.path.append('../../../../../../../inputFiles/hydraulicFracturing/scripts')
import hydrofractureFigure

hydrofractureFigure.main( geosDir='../../../../../../..',xmlFilePrefix='kgdToughnessDominated')
