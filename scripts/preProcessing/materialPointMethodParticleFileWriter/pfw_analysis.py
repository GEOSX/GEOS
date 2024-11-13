import numpy as np   
import pfw_geometryObjects as geom              
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import os
import re
import importlib
from cycler import cycler
import logging
from pfw_geometryObjects import countFileLines
import argparse


def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order"""
    with open(filename, 'rb') as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size)).decode(encoding='utf-8')
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first 
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment


# =================================================================
# DATA IMPORT/EXPORT
# =================================================================


def has_mpm_file(path):
  for p in os.listdir(path):
    if re.match(".*pfw_input_.*\.py", p):
      return True
  return False


def format_file_paths(runLocations, paths):
  files = []
  # Base condition path must be a directory, not a file
  for runLocation in runLocations:
    for pp in paths:
      p = runLocation + pp
      if os.path.isdir(p):
        if has_mpm_file(p):
          files.append(p)
        else:
          for x in formatFilePaths(runLocations, [p + "/" + y for y in os.listdir(p)]):
            files.append(x)

  return files


def read_from_reaction_file(filename):
  #First line are the headers so index starts at 1
  # time,F00,F11,F22,Rx-,Rx+,Ry-,Ry+,Rz-,Rz+
  react_data = np.genfromtxt(filename, delimiter=',')
  time = DataObj("Time", react_data[1:,0])
  F00 = DataObj("F00", react_data[1:,1])
  F11 = DataObj("F11", react_data[1:,2])
  F22 = DataObj("F22", react_data[1:,3])
  Lx = DataObj("Lx", react_data[1:,4])
  Ly = DataObj("Ly", react_data[1:,5])
  Lz = DataObj("Lz", react_data[1:,6])
  Rxm = DataObj("Rxm", react_data[1:,7])
  Rxp = DataObj("Rxp", react_data[1:,8])
  Rym = DataObj("Rym", react_data[1:,9])
  Ryp = DataObj("Ryp", react_data[1:,10])
  Rzm = DataObj("Rzm", react_data[1:,11])
  Rzp = DataObj("Rzp", react_data[1:,12])
  L00 = DataObj("L00", react_data[1:,13])
  L11 = DataObj("L11", react_data[1:,14])
  L22 = DataObj("L22", react_data[1:,15])

  return time, F00, F11, F22, Lx, Ly, Lz, Rxm, Rxp, Rym, Ryp, Rzm, Rzp


def read_from_box_average_file(filename):
  box_data = np.genfromtxt(filename, delimiter=',')

  time = DataObj("Time", box_data[1:,0])
  sxx = DataObj("Sxx", box_data[1:,1]) # xx
  syy = DataObj("Syy", box_data[1:,2]) # yy
  szz = DataObj("Szz", box_data[1:,3]) # zz 
  sxy = DataObj("Sxy", box_data[1:,4]) # xy
  syz = DataObj("Syz", box_data[1:,5]) # yz
  sxz = DataObj("Sxz", box_data[1:,6]) # xz
  density = DataObj("Density", box_data[1:,7]) # mass
  damage = DataObj("Dmg", box_data[1:,8]) # damage
  internalEnergy = DataObj("IE", box_data[1:,9]) # internal energy
  kineticEnergy = DataObj("KE", box_data[1:,10]) # kinetic energy
  epxx = DataObj("epxx", box_data[1:,11]) # plastic strain xx
  epyy = DataObj("epyy", box_data[1:,12]) # plastic strain xx
  epzz = DataObj("epzz", box_data[1:,13]) # plastic strain xx
  epyz = DataObj("epyz", box_data[1:,14]) # plastic strain xx
  epxz = DataObj("epxz", box_data[1:,15]) # plastic strain xx
  epxy = DataObj("epxy", box_data[1:,16]) # plastic strain xx
  matVol = DataObj("MatVol", box_data[1:,17]) # material volume

  return time, sxx, syy, szz, sxy, syz, sxz, density, damage, internalEnergy, kineticEnergy, epxx, epyy, epzz, epyz, epxz, epxy, matVol


def write_data_to_csv(filename, data_array):
  num_fields = len(data_array)
  headers = []
  num_entries = 0
  for n in range(num_fields):
    field = data_array[n]
    if n == 0:
        num_entries = len(field.getData())
    headers.append(field.name)
    assert num_entries == len(field.getData()), "Field " + field.name +  " had different numbe of entries. Expected of " + str(num_entries) + ". Got " + str(len(field.getData()))
        
  delimiter = ","
  
  with open(filename, 'w') as f:
    for h, header in enumerate(headers):
        f.write(header)
        if h != num_fields-1:
            f.write(delimiter)
    f.write("\n")
    for i in range(num_entries):
        for h in range(num_fields):
            f.write(("{:"+ data_array[h].format +"}").format(data_array[h].getData()[i]))
            if h != num_fields-1:
                f.write(delimiter)
        f.write("\n")


def write_data_to_console(filename, data_array):
  num_fields = len(data_array)
  headers = []
  num_entries = 0
  for n in range(num_fields):
    field = data_array[n]
    if n == 0:
      num_entries = len(field.getData())
    headers.append(field.name)
    assert num_entries == len(field.getData()), "Field " + field.name +  " had different numbe of entries. Expected of " + str(num_entries) + ". Got " + str(len(field.getData()))
        
  delimiter = ","
  
  for h, header in enumerate(headers):
    print(header, end="," if h != num_fields-1 else "\n")
  
  for i in range(num_entries):
    for h in range(num_fields):
      print(("{:"+ data_array[h].format +"}").format(data_array[h].getData()[i]),end="," if h != num_fields-1 else "\n")


# =================================================================
# END DATA IMPORT/EXPORT
# =================================================================


# =================================================================
# DATA POSTPROCESSING
# =================================================================


class Trim:
  def __init__(self, length):
    self.length = length
  
  def postprocess(self,x):
    return x[:self.length]


class RemoveNonMonotonicEntries:
  def __init__(self, x_in):
    self.x_in = x_in
    self.maxX = 0.0
    self.mask = np.ones(len(self.x_in), dtype=bool)
    for ii,t in enumerate(self.x_in):
      if (t<=self.maxX):
        self.mask[ii] = False
      else:
        self.maxX = t 

  def postprocess(self, x):
    return x[self.mask,...]


class MedianFilter:
  def __init__(self, window_size):
    self.window_size = window_size

  # Could probably use numpy matrices to speed this up
  def postprocess(self, x):
    x = np.array(x)
    N = len(x)
    x_out = np.copy(x)
    for n in range(self.window_size+1, N-self.window_size):
      x_out[n] = np.median(x[(n-self.window_size):(n+self.window_size+1)])
    return x_out 


class SubSample:
  def __init__(self, method="average", stride=None, numSamples=None):
    self.method = method
    self.stride = stride
    self.numSamples = numSamples

  def postprocess(self, x):
    if len(x) < self.numSamples:
      return x
    
    numX = len(x)
    if self.stride == None:
      self.stride = int(np.floor( numX / self.numSamples ))

    if self.method == "average":
      return np.average(np.pad(x, (0, self.stride - numX % self.stride), mode='constant', constant_values=x[-1]).reshape(-1, int(self.stride)), axis=1)

    if self.method == "nearest":
      return x[np.round(np.linspace(0, numX-1, num=self.numSamples)).astype(int)]

    # if method == "moving average":
    #   ret = np.cumsum(a, dtype=float)
    #   ret[stride:] = ret[stride:] - ret[:-stride]
    #   return ret[stride - 1:] / stride

    # match method:
    #   case "average":
    #     return np.average(np.pad(x, (0, numX % stride), mode='constant', constant_values=x[-1]).reshape(-1, stride), axis=1) 
    #   case _:
    #     print("No matching subsample method found")
    #     sys.exit(-1)
    

def compute_domain_strain(F00, F11, F22, engineeringStrain=False):
  if engineeringStrain:
      exx=F00-1.0
      eyy=F11-1.0
      ezz=F22-1.0
  else:
      exx=np.log(F00)
      eyy=np.log(F11)
      ezz=np.log(F22)

  J = F00*F11*F22

  return DataObj("exx", exx), DataObj("eyy", eyy), DataObj("ezz", ezz), DataObj("J", J)


def compute_domain_stress(Axx0, Ayy0, Azz0, F00, F11, F22, Rxm, Rxp, Rym, Ryp, Rzm, Rzp, engineeringStress=False):
  Ax=Axx0
  Ay=Ayy0
  Az=Azz0
  if not engineeringStress:
      Ax=Ax*F11*F22
      Ay=Ay*F00*F22
      Az=Az*F00*F11
  
  rsxm = -Rxm/Ax
  rsxp = Rxp/Ax
  rsxx=0.5*(Rxp-Rxm)/Ax
  rsym = -Rym/Ax
  rsyp = Ryp/Ax
  rsyy=0.5*(Ryp-Rym)/Ay
  rszm = -Rzm/Ax
  rszp = Rzp/Ax
  rszz=0.5*(Rzp-Rzm)/Az

  return DataObj("Rsxm", rsxm), DataObj("Rsxp", rsxp), DataObj("Rsxm", rsxx), DataObj("Rsym", rsym), DataObj("Rsyp", rsyp), DataObj("Rsyy", rsyy), DataObj("Rszm", rszm), DataObj("Rszp", rszp), DataObj("Rszz", rszz)


def compute_pressure(bsxx, bsyy, bszz):
  return DataObj("Pressure", -(bsxx + bsyy + bszz)/3)


# =================================================================
# END DATA POSTPROCESSING
# =================================================================


# =================================================================
# DATA VISUALIZATION/PLOTTING
# =================================================================


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


# =================================================================
# END DATA VISUALIZATION/PLOTTING
# =================================================================

class AnalysisOptions:
    def __init__(self):
        self.plotReactions = True # If false will only output results to console window
        self.readBoxSums = True # Read the boxAverage csv file output by jbos
        self.livePlot = False # Leave plotting window open, otherwise just saves image

        self.writePostprocessFile = True # writes data to file

        self.useEngineeringStrain=True

        self.flipCompression = True # Flips stress for compression simulations
        self.filterStresses = True # Removes spike artifacts in reaction data above a stress threshold
        self.stressMaxThreshold = np.inf # GPa
        self.stressMinThreshold = -stressMaxThreshold

        self.enforceAxesLimits = False # Handy for manually specifying axes limits in plots
        self.yStressMax=5
        self.yStressMin=-5

        self.filterByMedian = True
        self.windowSize = 5


class MPMJob:
  def __init__(self, job_dir_path):
      self.job_dir_path = job_dir_path
      self.job_name = os.path.basename(os.path.normpath(self.job_dir_path))
      self.job_input_file = 'pfw_input_'+ self.job_name + '.py'
      self.fields = {}
      self.gather_job_metadata()
      
  def gather_job_metadata(self):
      # Gather job meta data
      assert os.path.isfile(os.path.join(self.job_dir_path, self.job_input_file)), "Could not find " + self.job_input_file + " inside " + self.job_dir_path + ". Aborting..."
      
      sys.path.append(self.job_dir_path+"/")
      job = importlib.import_module(self.job_input_file[:-3])
      pfw = job.pfw

      # # Read stop time from XML
      # stop_time = pfw["endTime"] # This should be read from XML instead, in case someone modified the input or xml manually
      
      # # Read time from slurm outputs if available
      # lastTimestep = 0.0
      # slurmFile = "slurm-"+str(jobID)+".out"
      # for line in reverse_readline(slurmFile):
      #   if 'Time:' in line:
      #     line_terms = re.split(": |, |\s", line)
      #     # print(line_terms)
      #     lastTimestep = float(line_terms[1])
      #     # print("Last timestep was", lastTimestep)
      #     break
      # self.percent_complete = lastTimestep/stopTime

      self.density = job.density if hasattr( job, 'density') else -1 # Need to a default or message if no density is defined
      self.hasReactionFile = os.path.isfile(os.path.join(self.job_dir_path, "reactionHistory.csv"))
      self.hasBoxAverageFile = os.path.isfile(os.path.join(self.job_dir_path, "boxAverageHistory.csv"))

      self.domainWidth0 = job.domainWidth
      self.domainLength0 = job.domainLength
      self.domainHeight0 = job.domainHeight
      self.domainVolume0 = self.domainWidth0*self.domainLength0*self.domainHeight0

      self.sampleWidth = job.sampleWidth if hasattr( job, 'sampleWidth' ) else self.domainWidth0
      self.sampleHeight = job.sampleHeight if hasattr( job, 'sampleHeight' ) else self.domainHeight0
      self.sampleLength = job.sampleLength if hasattr( job, 'sampleLength' ) else self.domainLength0
      
      self.periodic = pfw["periodic"] if "periodic" in pfw else [False, False, False]

      self.NI = pfw["nI"]   # total grid cells in the x-direction (counting ghosts)
      self.NJ = pfw["nJ"]   # total grid cells in the y-direction (counting ghosts)
      self.NK = pfw["nK"]   # total grid cells in the z-direction (counting ghosts)
      self.ppc = pfw["ppc"] if hasattr(job, 'ppc') else 2 # particles per cell in each direction
      self.ppcx = pfw["ppcx"] if hasattr(job, 'ppcx') else self.ppc
      self.ppcy = pfw["ppcy"] if hasattr(job, 'ppcy') else self.ppc
      self.ppcz = pfw["ppcz"] if hasattr(job, 'ppcz') else self.ppc
      
      # interior discretizatiom
      self.nI = self.NI if self.periodic[0] else self.NI-2   # interior grid cells in the x-direction
      self.nJ = self.NJ if self.periodic[1] else self.NJ-2   # interior grid cells in the y-direction
      self.nK = self.NK if self.periodic[2] else self.NK-2   # interior grid cells in the z-direction

      self.xmin = pfw["xmin"]
      self.xmax = pfw["xmax"]
      self.ymin = pfw["ymin"]
      self.ymax = pfw["ymax"]
      self.zmin = pfw["zmin"]
      self.zmax = pfw["zmax"]

      # grid cell spacing
      self.dX = (self.xmax - self.xmin)/(self.nI * self.ppcx)
      self.dY = (self.ymax - self.ymin)/(self.nJ * self.ppcy)
      self.dZ = (self.zmax - self.zmin)/(self.nK * self.ppcz)
      
      # Eventually this should be generalized as an input for different deformation directions
      self.Axx0 = self.sampleHeight*self.sampleLength
      self.Ayy0 = self.sampleWidth*self.sampleLength
      self.Azz0 = self.sampleWidth*self.sampleHeight
      self.V0 = self.sampleWidth*self.sampleHeight*self.sampleLength

      self.numParticles = geom.countFileLines( os.path.join(self.job_dir_path, "mpmParticleFile_" + self.job_name) )

      self.planeStrain = job.planeStrain == 1 if hasattr( job, 'planeStrain' ) else False

      # Report metadata
      print("Meta data read from " + self.job_input_file + ":")
      # print("\tPercent Complete =", 100*self.percent_complete)
      print("\tDensity =", self.density)
      print("\thasReactionFile? =", self.hasReactionFile)
      print("\thasBoxAverageFile? =", self.hasBoxAverageFile)
      print("\tInitial Domain (Width, Height, Length) = (",self.domainWidth0,",", self.domainHeight0,",", self.domainLength0,")")
      print("\tPeriodic? = ", self.periodic)
      print("\tPlaneStrain? =", self.planeStrain)
      print("\tNum Particles =", self.numParticles)     

  def read_reaction_file(self):
    assert self.hasReactionFile, "No reaction file present!"

    #First line are the headers so index starts at 1
    # time,F00,F11,F22,Rx-,Rx+,Ry-,Ry+,Rz-,Rz+
    react_data = np.genfromtxt(os.path.join(self.job_dir_path, "reactionHistory.csv"), delimiter=',')
    self.registerField(DataObj("Time", react_data[1:,0]))
    self.registerField(DataObj("F00", react_data[1:,1]))
    self.registerField(DataObj("F11", react_data[1:,2]))
    self.registerField(DataObj("F22", react_data[1:,3]))
    self.registerField(DataObj("Lx", react_data[1:,4]))
    self.registerField(DataObj("Ly", react_data[1:,5]))
    self.registerField(DataObj("Lz", react_data[1:,6]))
    self.registerField(DataObj("Rxm", react_data[1:,7]))
    self.registerField(DataObj("Rxp", react_data[1:,8]))
    self.registerField(DataObj("Rym", react_data[1:,9]))
    self.registerField(DataObj("Ryp", react_data[1:,10]))
    self.registerField(DataObj("Rzm", react_data[1:,11]))
    self.registerField(DataObj("Rzp", react_data[1:,12]))
    self.registerField(DataObj("L00", react_data[1:,13]))
    self.registerField(DataObj("L11", react_data[1:,14]))
    self.registerField(DataObj("L22", react_data[1:,15]))

  def read_from_box_average_file(self):
    assert self.hasBoxAverageFile, "No box average file present!"

    box_data = np.genfromtxt(os.path.join(self.job_dir_path, "boxAverageHistory.csv"), delimiter=',')
    self.registerField(DataObj("BTime", box_data[1:,0]))
    self.registerField(DataObj("BSxx", box_data[1:,1])) # xx
    self.registerField(DataObj("BSyy", box_data[1:,2])) # yy
    self.registerField(DataObj("BSzz", box_data[1:,3])) # zz 
    self.registerField(DataObj("BSxy", box_data[1:,4])) # xy
    self.registerField(DataObj("BSyz", box_data[1:,5])) # yz
    self.registerField(DataObj("BSxz", box_data[1:,6])) # xz
    self.registerField(DataObj("BDensity", box_data[1:,7])) # mass
    self.registerField(DataObj("BDmg", box_data[1:,8])) # damage
    self.registerField(DataObj("BIE", box_data[1:,9])) # internal energy
    self.registerField(DataObj("BKE", box_data[1:,10])) # kinetic energy
    self.registerField(DataObj("Bepxx", box_data[1:,11])) # plastic strain xx
    self.registerField(DataObj("Bepyy", box_data[1:,12])) # plastic strain xx
    self.registerField(DataObj("Bepzz", box_data[1:,13])) # plastic strain xx
    self.registerField(DataObj("Bepyz", box_data[1:,14])) # plastic strain xx
    self.registerField(DataObj("Bepxz", box_data[1:,15])) # plastic strain xx
    self.registerField(DataObj("Bepxy", box_data[1:,16])) # plastic strain xx
    self.registerField(DataObj("BMatVol", box_data[1:,17])) # material volume

  def registerField(self, field):
    self.fields[field.name] = field

  def applyPostProcess(self, fieldname, filter):
    if fieldname == "all":
      for name, field in self.fields.items():
        field.applyPostProcess(filter)
    else:
      self.fields[fieldname].applyPostProcess(filter)


  def compute_domain_strain(self, engineeringStrain=False):
    if engineeringStrain:
        exx=self.fields["F00"].getData()-1.0
        eyy=self.fields["F11"].getData()-1.0
        ezz=self.fields["F22"].getData()-1.0
    else:
        exx=np.log(self.fields["F00"].getData())
        eyy=np.log(self.fields["F11"].getData())
        ezz=np.log(self.fields["F22"].getData())

    J =self.fields["F00"].getData() *  self.fields["F11"].getData() * self.fields["F22"].getData()

    self.registerField(DataObj("exx", exx))
    self.registerField(DataObj("eyy", eyy))
    self.registerField(DataObj("ezz", ezz))
    self.registerField(DataObj("J", J))


  def compute_domain_stress(self, engineeringStress=False):
    Ax=self.Axx0
    Ay=self.Ayy0
    Az=self.Azz0
    if not engineeringStress:
        Ax=Ax*self.fields["F11"].getData()*self.fields["F22"].getData()
        Ay=Ay*self.fields["F00"].getData()*self.fields["F22"].getData()
        Az=Az*self.fields["F11"].getData()*self.fields["F00"].getData()
    
    self.registerField(DataObj("Rsxm", -self.fields["Rxm"].getData()/Ax))
    self.registerField(DataObj("Rsxp", self.fields["Rxp"].getData()/Ax))
    self.registerField(DataObj("Rsxx", 0.5*(self.fields["Rxp"].getData()-self.fields["Rxm"].getData())/Ax))
    self.registerField(DataObj("Rsym", -self.fields["Rym"].getData()/Ax))
    self.registerField(DataObj("Rsyp", self.fields["Ryp"].getData()/Ax))
    self.registerField(DataObj("Rsyy", 0.5*(self.fields["Ryp"].getData()-self.fields["Rym"].getData())/Ay))
    self.registerField(DataObj("Rszm", -self.fields["Rzm"].getData()/Ax))
    self.registerField(DataObj("Rszp", self.fields["Rzp"].getData()/Ax))
    self.registerField(DataObj("Rszz", 0.5*(self.fields["Rzp"].getData()-self.fields["Rzm"].getData())/Az))   
        

class DataObj:
    def __init__(self, name, data, format=".6f"):
      self.name = name
      self.data = data
      self.format =format
      self.processedData = self.data

    def applyPostProcess(self, filter):
      self.processedData = filter.postprocess(self.processedData)
    
    def getData(self):
      return self.processedData


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Analyze GEOS MPM Jobs')
  parser.add_argument('jobdir', help="location of job to analyze")
  parser.add_argument('-i', '--interactive', action='store_true', help="display plot interactively")
  parser.add_argument('-c','--console',action='store_true', default=False, help="flag to write output to console")
  parser.add_argument('-x', '--xyz', default="xyz", help="list of directions to plot reactions and displacements")
  parser.add_argument('-e', '--export', action='store_true', default=False, help="write output to csv")
  parser.add_argument('-s', '--save', default=False, help="flag to save plot to png file")
  parser.add_argument('-p', '--plot', default=None, help="plot job data (Default reactions and FTable)")
  parser.add_argument('-f','--fields', nargs='+', help="list of field names to output", type=str)
  parser.add_argument('-o', '--output', default="out", help="name of output files")
  args = parser.parse_args()
  print("Command line args",args)

  assert os.path.isdir(args.jobdir), "Could not find job directory with name:" + args.jobdir

  job = MPMJob(args.jobdir)

  if args.plot is not None:
    if args.plot=="reactions":
      job.read_reaction_file()

      fig = plt.figure()
      ax = fig.add_subplot()
      ax.plot(job.fields["Time"].getData(), job.fields["Ryp"].getData())
      plt.show()
