# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 09:00:00 2017
Updated on Wed Apr 17 10:30:00 2024
@author: homel1, crook5
Geometry object functions for the particle file writer.
"""
import numpy as np                   
from sklearn.neighbors import KDTree         
from random import seed
from random import random
import math
from scipy.optimize import minimize_scalar
from scipy.spatial import Voronoi, ConvexHull
from abc import ABCMeta, abstractmethod
import scipy
np.random.seed(1)

# ===========================================
# DEFAULTS
# ===========================================

_defaultMat = 0
_defaultGroup = 0
_defaultVelocity = np.array([0.0, 0.0, 0.0])
_defaultDamage = 0.0
_defaultPorosity = 0.0
_defaultTemperature = 300.0
_defaultSurfaceFlag = 2 # Basic surface flag (0 = internal, 1=fully damaged, 2=normal surface flag, 3=cohesive)
_defaultSurfaceNormal = np.array([0.0, 0.0, 0.0]) # 0 vector turns off surface normal use in solver (defaults to implicit surface normals)
_defaultSurfacePosition = np.array([0.0, 0.0, 0.0])
_defaultMatDir = np.array([1.0, 0.0, 0.0])
_defaultSurfaceTraction = np.array([0.0,0.0,0.0])
_defaultParticleType = 2 # 0 (Single point), 1 (Single Point Bsplines), 2 (CPDI), 3 (CPTI), 4 (CPDI2)

# END DEFAULTS
# ===========================================


# ===========================================
_gyroid_rhoVsC = np.array([[0.0000000000000000, 0.0000000000000000],
                          [0.0097600000000000, 0.0158666295635848],
                          [0.0204160000000000, 0.0317332591271696],
                          [0.0304960000000000, 0.0475998886907544],
                          [0.0399040000000000, 0.0634665182543393],
                          [0.0499840000000000, 0.0793331478179241],
                          [0.0604000000000000, 0.0951997773815089],
                          [0.0710240000000000, 0.1110664069450940],
                          [0.0810080000000000, 0.1269330365086790],
                          [0.0908480000000000, 0.1427996660722630],
                          [0.1011200000000000, 0.1586662956358480],
                          [0.1116800000000000, 0.1745329251994330],
                          [0.1209920000000000, 0.1903995547630180],
                          [0.1322400000000000, 0.2062661843266030],
                          [0.1411680000000000, 0.2221328138901870],
                          [0.1519200000000000, 0.2379994434537720],
                          [0.1621920000000000, 0.2538660730173570],
                          [0.1721760000000000, 0.2697327025809420],
                          [0.1821120000000000, 0.2855993321445270],
                          [0.1924480000000000, 0.3014659617081110],
                          [0.2032960000000000, 0.3173325912716960],
                          [0.2130880000000000, 0.3331992208352810],
                          [0.2232160000000000, 0.3490658503988660],
                          [0.2331040000000000, 0.3649324799624510],
                          [0.2442080000000000, 0.3807991095260360],
                          [0.2541920000000000, 0.3966657390896200],
                          [0.2642720000000000, 0.4125323686532050],
                          [0.2740640000000000, 0.4283989982167900],
                          [0.2852000000000000, 0.4442656277803750],
                          [0.2943680000000000, 0.4601322573439600],
                          [0.3052800000000000, 0.4759988869075440],
                          [0.3150240000000000, 0.4918655164711290],
                          [0.3260160000000000, 0.5077321460347140],
                          [0.3362400000000000, 0.5235987755982990],
                          [0.3458880000000000, 0.5394654051618840],
                          [0.3569280000000000, 0.5553320347254680],
                          [0.3675520000000000, 0.5711986642890530],
                          [0.3776800000000000, 0.5870652938526380],
                          [0.3878560000000000, 0.6029319234162230],
                          [0.3992320000000000, 0.6187985529798080],
                          [0.4091680000000000, 0.6346651825433930],
                          [0.4196800000000000, 0.6505318121069770],
                          [0.4294880000000000, 0.6663984416705620],
                          [0.4405280000000000, 0.6822650712341470],
                          [0.4508480000000000, 0.6981317007977320],
                          [0.4613600000000000, 0.7139983303613170],
                          [0.4727360000000000, 0.7298649599249010],
                          [0.4826400000000000, 0.7457315894884860],
                          [0.4941600000000000, 0.7615982190520710],
                          [0.5037600000000000, 0.7774648486156560],
                          [0.5157120000000000, 0.7933314781792410],
                          [0.5259360000000000, 0.8091981077428260],
                          [0.5365120000000000, 0.8250647373064100],
                          [0.5471200000000000, 0.8409313668699950],
                          [0.5584480000000000, 0.8567979964335800],
                          [0.5691040000000000, 0.8726646259971650],
                          [0.5805760000000000, 0.8885312555607500],
                          [0.5908640000000000, 0.9043978851243340],
                          [0.6026240000000000, 0.9202645146879190],
                          [0.6137600000000000, 0.9361311442515040],
                          [0.6247040000000000, 0.9519977738150890],
                          [0.6359360000000000, 0.9678644033786740],
                          [0.6479040000000000, 0.9837310329422580],
                          [0.6594720000000000, 0.9995976625058430],
                          [0.6711360000000000, 1.0154642920694300],
                          [0.6830400000000000, 1.0313309216330100],
                          [0.6943360000000000, 1.0471975511966000],
                          [0.7055680000000000, 1.0630641807601800],
                          [0.7171840000000000, 1.0789308103237700],
                          [0.7285120000000000, 1.0947974398873500],
                          [0.7403840000000000, 1.1106640694509400],
                          [0.7519040000000000, 1.1265306990145200],
                          [0.7631360000000000, 1.1423973285781100],
                          [0.7747520000000000, 1.1582639581416900],
                          [0.7862400000000000, 1.1741305877052800],
                          [0.7979520000000000, 1.1899972172688600],
                          [0.8096160000000000, 1.2058638468324500],
                          [0.8206560000000000, 1.2217304763960300],
                          [0.8325280000000000, 1.2375971059596200],
                          [0.8447680000000000, 1.2534637355232000],
                          [0.8568640000000000, 1.2693303650867900],
                          [0.8688320000000000, 1.2851969946503700],
                          [0.8807840000000000, 1.3010636242139500],
                          [0.8925920000000000, 1.3169302537775400],
                          [0.9049440000000000, 1.3327968833411200],
                          [0.9166560000000000, 1.3486635129047100],
                          [0.9291840000000000, 1.3645301424682900],
                          [0.9419200000000000, 1.3803967720318800],
                          [0.9544000000000000, 1.3962634015954600],
                          [0.9667520000000000, 1.4121300311590500],
                          [0.9769760000000000, 1.4279966607226300],
                          [0.9847200000000000, 1.4438632902862200],
                          [0.9913920000000000, 1.4597299198498000],
                          [0.9959680000000000, 1.4755965494133900],
                          [0.9992000000000000, 1.4914631789769700],
                          [1.0000000000000000, 1.5073298085405600]])

_schwarzDiamond_rhoVsC = np.array([[0.000000000000000, 0.0000000000000000],
                                  [0.019200000000000, 0.0158666295635848],
                                  [0.040416000000000, 0.0317332591271696],
                                  [0.057728000000000, 0.0475998886907544],
                                  [0.076832000000000, 0.0634665182543393],
                                  [0.095936000000000, 0.0793331478179241],
                                  [0.114272000000000, 0.0951997773815089],
                                  [0.133504000000000, 0.1110664069450940],
                                  [0.151264000000000, 0.1269330365086790],
                                  [0.169600000000000, 0.1427996660722630],
                                  [0.189088000000000, 0.1586662956358480],
                                  [0.208128000000000, 0.1745329251994330],
                                  [0.225792000000000, 0.1903995547630180],
                                  [0.244800000000000, 0.2062661843266030],
                                  [0.263328000000000, 0.2221328138901870],
                                  [0.281024000000000, 0.2379994434537720],
                                  [0.301184000000000, 0.2538660730173570],
                                  [0.317600000000000, 0.2697327025809420],
                                  [0.337376000000000, 0.2855993321445270],
                                  [0.355936000000000, 0.3014659617081110],
                                  [0.373216000000000, 0.3173325912716960],
                                  [0.391840000000000, 0.3331992208352810],
                                  [0.411136000000000, 0.3490658503988660],
                                  [0.429696000000000, 0.3649324799624510],
                                  [0.446880000000000, 0.3807991095260360],
                                  [0.465024000000000, 0.3966657390896200],
                                  [0.485088000000000, 0.4125323686532050],
                                  [0.501824000000000, 0.4283989982167900],
                                  [0.522368000000000, 0.4442656277803750],
                                  [0.540704000000000, 0.4601322573439600],
                                  [0.559232000000000, 0.4759988869075440],
                                  [0.577888000000000, 0.4918655164711290],
                                  [0.596512000000000, 0.5077321460347140],
                                  [0.616672000000000, 0.5235987755982990],
                                  [0.634912000000000, 0.5394654051618840],
                                  [0.651552000000000, 0.5553320347254680],
                                  [0.671808000000000, 0.5711986642890530],
                                  [0.690528000000000, 0.5870652938526380],
                                  [0.710144000000000, 0.6029319234162230],
                                  [0.728672000000000, 0.6187985529798080],
                                  [0.748064000000000, 0.6346651825433930],
                                  [0.768736000000000, 0.6505318121069770],
                                  [0.786112000000000, 0.6663984416705620],
                                  [0.805600000000000, 0.6822650712341470],
                                  [0.825088000000000, 0.6981317007977320],
                                  [0.843840000000000, 0.7139983303613170],
                                  [0.859456000000000, 0.7298649599249010],
                                  [0.874016000000000, 0.7457315894884860],
                                  [0.888448000000000, 0.7615982190520710],
                                  [0.900480000000000, 0.7774648486156560],
                                  [0.911104000000000, 0.7933314781792410],
                                  [0.922560000000000, 0.8091981077428260],
                                  [0.932800000000000, 0.8250647373064100],
                                  [0.941952000000000, 0.8409313668699950],
                                  [0.950720000000000, 0.8567979964335800],
                                  [0.959712000000000, 0.8726646259971650],
                                  [0.966944000000000, 0.8885312555607500],
                                  [0.973728000000000, 0.9043978851243340],
                                  [0.979616000000000, 0.9202645146879190],
                                  [0.985120000000000, 0.9361311442515040],
                                  [0.990368000000000, 0.9519977738150890],
                                  [0.994400000000000, 0.9678644033786740],
                                  [0.998016000000000, 0.9837310329422580],
                                  [0.999936000000000, 0.9995976625058430],
                                  [1.000000000000000, 1.0154642920694300]])

_schwarzPrimitive_rhoVsC = np.array([[0.000000000000000, 0.0000000000000000],
                                    [0.017808000000000, 0.0317332591271696],
                                    [0.036560000000000, 0.0634665182543393],
                                    [0.053792000000000, 0.0951997773815089],
                                    [0.072536000000000, 0.1269330365086790],
                                    [0.089680000000000, 0.1586662956358480],
                                    [0.108160000000000, 0.1903995547630180],
                                    [0.125728000000000, 0.2221328138901870],
                                    [0.143952000000000, 0.2538660730173570],
                                    [0.161592000000000, 0.2855993321445270],
                                    [0.180024000000000, 0.3173325912716960],
                                    [0.197624000000000, 0.3490658503988660],
                                    [0.215936000000000, 0.3807991095260360],
                                    [0.234152000000000, 0.4125323686532050],
                                    [0.251992000000000, 0.4442656277803750],
                                    [0.270160000000000, 0.4759988869075440],
                                    [0.287896000000000, 0.5077321460347140],
                                    [0.306312000000000, 0.5394654051618840],
                                    [0.324120000000000, 0.5711986642890530],
                                    [0.342648000000000, 0.6029319234162230],
                                    [0.360176000000000, 0.6346651825433930],
                                    [0.378824000000000, 0.6663984416705620],
                                    [0.396680000000000, 0.6981317007977320],
                                    [0.415192000000000, 0.7298649599249010],
                                    [0.433024000000000, 0.7615982190520710],
                                    [0.451384000000000, 0.7933314781792410],
                                    [0.469320000000000, 0.8250647373064100],
                                    [0.488088000000000, 0.8567979964335800],
                                    [0.505992000000000, 0.8885312555607500],
                                    [0.525176000000000, 0.9202645146879190],
                                    [0.542528000000000, 0.9519977738150890],
                                    [0.562360000000000, 0.9837310329422580],
                                    [0.579016000000000, 1.0154642920694300],
                                    [0.595840000000000, 1.0471975511966000],
                                    [0.610416000000000, 1.0789308103237700],
                                    [0.625440000000000, 1.1106640694509400],
                                    [0.638520000000000, 1.1423973285781100],
                                    [0.651656000000000, 1.1741305877052800],
                                    [0.664400000000000, 1.2058638468324500],
                                    [0.676448000000000, 1.2375971059596200],
                                    [0.687976000000000, 1.2693303650867900],
                                    [0.699256000000000, 1.3010636242139500],
                                    [0.710736000000000, 1.3327968833411200],
                                    [0.721248000000000, 1.3645301424682900],
                                    [0.731832000000000, 1.3962634015954600],
                                    [0.741632000000000, 1.4279966607226300],
                                    [0.751568000000000, 1.4597299198498000],
                                    [0.760976000000000, 1.4914631789769700],
                                    [0.770440000000000, 1.5231964381041400],
                                    [0.779320000000000, 1.5549296972313100],
                                    [0.788280000000000, 1.5866629563584800],
                                    [0.796848000000000, 1.6183962154856500],
                                    [0.805152000000000, 1.6501294746128200],
                                    [0.813176000000000, 1.6818627337399900],
                                    [0.821600000000000, 1.7135959928671600],
                                    [0.829120000000000, 1.7453292519943300],
                                    [0.836536000000000, 1.7770625111215000],
                                    [0.843928000000000, 1.8087957702486700],
                                    [0.851016000000000, 1.8405290293758400],
                                    [0.857976000000000, 1.8722622885030100],
                                    [0.864560000000000, 1.9039955476301800],
                                    [0.871400000000000, 1.9357288067573500],
                                    [0.877672000000000, 1.9674620658845200],
                                    [0.884200000000000, 1.9991953250116900],
                                    [0.890080000000000, 2.0309285841388600],
                                    [0.896184000000000, 2.0626618432660300],
                                    [0.901752000000000, 2.0943951023932000],
                                    [0.907424000000000, 2.1261283615203700],
                                    [0.912848000000000, 2.1578616206475300],
                                    [0.918208000000000, 2.1895948797747000],
                                    [0.923344000000000, 2.2213281389018700],
                                    [0.928224000000000, 2.2530613980290400],
                                    [0.932832000000000, 2.2847946571562100],
                                    [0.937784000000000, 2.3165279162833800],
                                    [0.942248000000000, 2.3482611754105500],
                                    [0.946624000000000, 2.3799944345377200],
                                    [0.950896000000000, 2.4117276936648900],
                                    [0.955056000000000, 2.4434609527920600],
                                    [0.958928000000000, 2.4751942119192300],
                                    [0.962624000000000, 2.5069274710464000],
                                    [0.966304000000000, 2.5386607301735700],
                                    [0.969856000000000, 2.5703939893007400],
                                    [0.973176000000000, 2.6021272484279100],
                                    [0.976304000000000, 2.6338605075550800],
                                    [0.979400000000000, 2.6655937666822500],
                                    [0.982144000000000, 2.6973270258094200],
                                    [0.984960000000000, 2.7290602849365900],
                                    [0.987488000000000, 2.7607935440637600],
                                    [0.989872000000000, 2.7925268031909300],
                                    [0.992208000000000, 2.8242600623181000],
                                    [0.994184000000000, 2.8559933214452700],
                                    [0.995872000000000, 2.8877265805724400],
                                    [0.997464000000000, 2.9194598396996100],
                                    [0.998800000000000, 2.9511930988267800],
                                    [0.999728000000000, 2.9829263579539500],
                                    [1.000000000000000, 3.0146596170811100]])

# ===========================================

# ===========================================
# UTILITY FUNCTIONS
# ===========================================

log_file = "job_file_log"
def log2file(msg):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")

    msg = current_time + " Rank " + str(g_rank) + ": " + msg + "\n"
    mode = MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND

    file_handle = MPI.File.Open(MPI.COMM_SELF, log_file + "_" + str(g_rank), mode)
    file_handle.Set_atomicity(True)

    b = bytearray()
    b.extend(map(ord, msg))

    file_handle.Write(b)
    file_handle.Sync()
    file_handle.Close()


def print2file(file_name, text):
    with open(file_name, 'a') as f:
        f.write(text + "\n")


def countFileLines(fname):
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b

    with open(fname, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count

def fOffsetLineNum(fname, line_number):
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b

    if line_number == 0:
        return 0

    count = 0
    with open(fname, "rb") as f:
        for buf in _make_gen(f.raw.read):
            count += buf.count(b"\n")
            if count >= line_number:
                chunk_pos = f.tell()
                break

        # Back track to line_number
        for ci in range(len(buf)-1, -1, -1):
            c = buf[ci]
            # print("Character ",c, )
            chunk_pos -= 1
            if c == 10: #b"\n":
                if count == line_number:
                    return chunk_pos + 1 # Add one char because first char of new line immediately follows new line char
                count -= 1


#Might be a faster way to do this
def fileOffsetFromLineNum(fname, line_number):
    char_pos = 0
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            if i == line_number:
                offset = len(line)
                return char_pos
            else:
                char_pos += len(line)


# used for periodicity
def mapToRange(x, xmin, xmax):
    # Make sure all inputs are np arrays
    xmin = np.array(xmin)
    xmax = np.array(xmax)
    dx = xmax-xmin
    x = np.array(x) - xmin

    x = x * (x < dx) + (x >= dx) * (np.fmod(x, dx))
    x = x * (x >= 0) + (x < 0) * (dx - np.fmod(np.absolute(x), dx))

    return x + xmin


# Smooth sigmoid interpolation between rf and rf0
def smoothStep3(r, rf, rf0):
  if r >= rf0:
    return 0

  if r <= rf:
    return 1

  r = (rf0-r)/(rf0-rf)

  return r*r*(3-2*r)


def inside_box(x, x0=np.array([0.0,0.0,0.0]), dx=np.array([1.0,1.0,1.0]), periodic=[False, False, False]):
  x = np.array(x)
  return np.all(np.logical_or(periodic, np.logical_and(x >= x0, x < x0+dx)))


# Checks if point x lies outside of box with minimum corner at x0 and dimensions dx
def outside(x, x0=np.array([0.0,0.0,0.0]), dx=np.array([1.0,1.0,1.0]), periodic=[False, False, False]):
    x = periodic * (np.fmod(x - x0, dx) + x0) + np.logical_not(periodic) * x
    return np.any(np.logical_or(x < x0, x > x0 + dx)), x
    

# Grid object used for efficiently performing random sequential absorption and poisson disc sampling
class Grid:
    def __init__(self, cell_size, **kwargs):
        self.cell_size = cell_size
        self.dim = kwargs.get('dim', 3)
        self.x0 = np.array(kwargs.get('x0', np.zeros(self.dim)))
        self.dx = np.array(kwargs.get('dx', np.ones(self.dim)))
        self.periodic = kwargs.get('periodic', [True for k in range(self.dim)])
        self.num_cells = np.ceil(self.dx / self.cell_size).astype(int)
        self.has_fractional_cells = np.multiply(self.num_cells, self.dx) > self.dx
        self.cells = -np.ones(tuple(self.num_cells)).astype(int)  # Cell value of -1 indicates empty cell

    def outside(self, x):
        return outside(x, self.x0, self.dx, self.periodic)

    def pos2cell(self, x):
        return np.floor((x - self.x0) / self.cell_size).astype(int)

    def place(self, cell, x):
        self.cells[tuple(cell)] = x
        return self

    def at(self, cell):
        return self.cells[tuple(cell)]

    def occupied(self, cell):
        return self.cells[tuple(cell)] >= 0

    def cells_in_range(self, cell, window):
        search_min = cell - window
        search_max = cell + window

        # If direction doesn't have periodic boundaries then cap search window to min and max cells
        # If grid has fractional cells and candidate cell is within range of them add 1 to search direction
        search_min = (search_min - np.ones(self.dim).astype(int) * np.logical_and(self.has_fractional_cells, cell < 1)) * self.periodic + np.logical_not(self.periodic) * np.maximum(
            np.zeros(self.dim).astype(int), search_min)

        search_max = (search_max + np.ones(self.dim).astype(int) * np.logical_and(self.has_fractional_cells, cell >= self.num_cells - 3)) * self.periodic + np.logical_not(self.periodic) * np.minimum(
            self.num_cells-1, search_max)

        search_range = search_max-search_min
        larger_than_cells = search_range > ( self.num_cells - 1 )

        search_min = np.zeros(self.dim).astype(int) * larger_than_cells + np.logical_not(larger_than_cells) * search_min
        search_max = (self.num_cells - 1) * larger_than_cells + np.logical_not(larger_than_cells) * search_max

        # Create list of all cells indices in search area
        cells_to_search_ranges = [np.array([]) for k in range(self.dim)]
        for k in range(self.dim):
            cells_to_search_ranges[k] = np.array(list(range(search_min[k], search_max[k] + 1))).astype(int)
        cells_to_search = combvec(cells_to_search_ranges).astype(int)


        cell_vals = []
        for c in cells_to_search:
            # Check if c index needs to be wrapped in periodic direction
            c = wrap(c, np.zeros(self.dim), self.num_cells - 1)

            val = self.cells[tuple(c)]
            if val >= 0:
                cell_vals.append(val)

        return cell_vals


# wraps integer index i to range of min and max (both inclusive)
# accepts nd arrays
def wrap(i, lower, upper):
    # Make sure all inputs are np arrays
    i = np.array(i).astype(int)
    lower = np.array(lower).astype(int)
    upper = np.array(upper).astype(int)

    r = upper - lower + 1
    i = i - lower

    i = i * (i < r) + (i >= r) * (lower + np.mod(i, r))
    i = i * (i >= 0) + (i < 0) * (upper - np.mod(np.absolute(i + 1), r))

    return i.astype(int)


# Generates all possible combinations of elements from list of offsets for each direction
def combvec(offsets):
    dim = len(offsets)

    combs = np.reshape(offsets[0], (len(offsets[0]), 1))
    for i in range(1, dim):
        new_combs = np.empty((0, i+1))
        for j in offsets[i]:
            for k in combs:
                new_combs = np.append(new_combs, np.reshape(np.append(k, j), (1, len(k)+1)), axis=0)

        combs = new_combs

    return combs


# Randomly samples points on surface of unit n-sphere
def random_direction(dim=3):
    direction = np.random.normal(size=dim)
    return direction / np.linalg.norm(direction)


# Adds images of pores for nearest neighbors
def add_images(pores, **kwargs):
    # Setup parameters
    dim = kwargs.get('dim', 3)
    x0 = np.array(kwargs.get('x0', np.zeros(dim)))
    dx = np.array(kwargs.get('dx', np.ones(dim)))
    periodic = kwargs.get('periodic', [True for k in range(dim)])

    num_pores = len(pores)
    boundary_pores = np.empty((0, dim + 1))
    offsets = [np.array([0, -1, 1]) for i in range(dim)]
    offset_combs = combvec(offsets)
    for ff in offset_combs:
        for p in pores:
            boundary_pores = np.append(boundary_pores, np.reshape(p + np.append(np.multiply(ff, dx), [0]), (1, dim + 1)), axis=0)

    return boundary_pores


# Adds pores images of those that lie on the boundary of the RVE
# x0 = corner with minimum x, y, and z coordinates
# dx = width, height, and length of domain
# pores = list of pores, format = x, y, z, r
def add_pore_images(pores, **kwargs):
  # Setup parameters
  dim = kwargs.get('dim', 3)
  x0 = np.array(kwargs.get('x0', np.zeros(dim)))
  dx = np.array(kwargs.get('dx', np.ones(dim)))
  periodic = kwargs.get('periodic', [True for k in range(dim)])

  num_pores = len(pores)

  indices = []
  boundary_pores = np.empty((0, dim + 1))
  for i in range(num_pores):
      offsets = [np.array([0]) for k in range(dim)]

      for k in range(dim):
          if periodic[k]:
              # Pore on - boundary
              if pores[i, k] - pores[i, dim] < x0[k]:
                  offsets[k] = np.append(offsets[k], 1.0)

              # Pore on + boundary
              if pores[i, k] + pores[i, dim] > x0[k] + dx[k]:
                  offsets[k] = np.append(offsets[k], -1.0)

      offset_combs = combvec(offsets)
      for ff in offset_combs:
          indices.append(i)
          boundary_pores = np.append(boundary_pores, np.reshape(pores[i, :] + np.append(np.multiply(ff, dx), [0]), (1, dim+1)), axis=0)

  return boundary_pores, indices


# Imports pores from a LAMMPS dump file
# pack_file = location of pack.dump file
# offset = translation of coordinates for simulation (default [0.0, 0.0, 0.0])
# add_boundary_images = boolean whether to include pores that lie on boundaries on other periodic surface 
# (e.g for generating the foam meshes in the particleFileWriter)
def import_pores(pack_file, **kwargs):
    # Setup parameters
    dim = kwargs.get('dim', 3)

    # x0 = np.array(kwargs.get('x0', np.zeros(dim)))
    # dx = np.array(kwargs.get('dx', np.ones(dim)))
    offset = np.array(kwargs.get('offset', np.zeros(dim)))
    scale = kwargs.get('scale', 1.0)
    add_boundary_images = kwargs.get('add_boundary_images', True)
    radius_scale = kwargs.get('radius_scale', 1.0)

    print("Reading pores from " + pack_file)

    file = open(pack_file, "r")
    lines = file.readlines()
    file.close()

    # Stores current data block being processed and list of voids
    curr_block = ""
    pores = []
    box = []
    curr_iter = 0
    num_atoms = 0
    count = 0
    for line in lines:
        count += 1
        # print("Line{}: {}".format(count, line.strip()))
        if line == "" or line == "\n":
            continue

        # Split string by spaces
        line_terms = line.split()

        # Determined what to do based on first entry in line (ITEM: = some header info, number
        if line_terms[0] == "ITEM:":
            curr_block = line_terms[1]

            if curr_block == "BOX":
                # For box get periodicity from same line
                periodic = [line_terms[3] == 'pp', line_terms[4] == 'pp', line_terms[5] == 'pp']
        elif curr_block == "ITERATION":
            curr_iter = int(line)
        elif curr_block == "NUMBER":
            num_atoms = int(line)
        elif curr_block == "BOX":
            box.append([float(i) for i in line_terms])
        elif curr_block == "ATOMS":
            atom = [float(x) for x in line_terms]
            # Output format radius, x, y, z
            pores.append([scale*(atom[2] + offset[0]), scale*(atom[3] + offset[1]), scale*(atom[4] + offset[2]), scale * radius_scale * atom[5]])

    pores = np.array(pores)  # convert python array to numpy array

    x0 = [scale*(box[0][0] + offset[0]), scale*(box[1][0] + offset[1]), scale*(box[2][0] + offset[2])]
    dx = [scale*box[0][1], scale*box[1][1], scale*box[2][1]]

    print("Box after offset and scaling:")
    print("\t" + str(x0[0]) + " " + str(dx[0]))
    print("\t" + str(x0[1]) + " " + str(dx[1]))
    print("\t" + str(x0[2]) + " " + str(dx[2]))
    print("Periodicity: " + str(periodic[0]) + " " + str(periodic[1]) + " " + str(periodic[2]))
    print("Pores from file: " + str(len(pores)))
    if add_boundary_images:
        output_pores = add_pore_images(pores, x0=x0, dx=dx, dim=dim, periodic=periodic)
        print("Pore with images: " + str(len(output_pores)))
    else:
        output_pores = pores

    # Reformt pores so radius is first, then x, y, z
    radii = output_pores[:,dim].copy()
    output_pores[:,1:dim+1] = output_pores[:,0:dim] # shift x, y, z right one column
    output_pores[:,0] = radii

    return np.array(output_pores)


# Generates densely packed random points in box domain with minimum separation
# Generalized for 2 and 3 dimensions
def poisson(spacing, **kwargs):
    # Setup parameters
    seed = kwargs.get('seed', np.random.randint(0, 1e8))
    dim = kwargs.get('dim', 3)
    x0 = np.array(kwargs.get('x0', np.zeros(dim)))
    dx = np.array(kwargs.get('dx', np.ones(dim)))
    trials = kwargs.get('trials', 30)
    periodic = kwargs.get('periodic', [True for k in range(dim)])

    x0 = x0[0:dim]
    dx = dx[0:dim]

    # Grid cell size if determined given longest cell size diagonal is radius of spacing
    grid = Grid(spacing / math.sqrt(dim), dim=dim, x0=x0, dx=dx, periodic=periodic)

    # First seed point is at center of domain
    seed_pts = np.array([dx/2 + x0])  # Check that this makes the right size array
    pores = np.empty((0, dim+1))

    # Set seed for random number generation
    np.random.seed(seed)

    # Terminates when there are no more seed points left
    while seed_pts.shape[0] != 0:
        # Get the first point in the seed array
        curr_seed_pt = seed_pts[0, :]
        seed_pts = np.delete(seed_pts, 0, axis=0)

        num_rejected = 0
        while num_rejected < trials:
            # Pick random distance and random angle from seed point
            d = spacing*(1 + np.random.random() * 0.1)  # random spacing between seed and candidate point
            direction = random_direction(dim)
            candidate_pt = curr_seed_pt + d * direction

            accept = True  # Accept candidate point unless found otherwise

            # If periodic boundaries are not enabled test point must be inside domain
            is_outside, candidate_pt = grid.outside(candidate_pt)
            if is_outside:
                num_rejected = num_rejected + 1
                continue

            candidate_cell = grid.pos2cell(candidate_pt)
            if grid.occupied(candidate_cell):
                accept = False
            else:
                neighbor_pores = grid.cells_in_range(candidate_cell, 2*np.ones(dim).astype(int))

                # print( "Neighbors", neighbor_pores)
                # for n in neighbor_pores:
                #     print(n, pores[n,:-1], grid.pos2cell(pores[n,:-1]))
                # print("Done Neighbors")

                # Check every cell for occupant, if any are within spacing reject candidate point
                for n in neighbor_pores:

                    # Check distance between candidate point and point already there
                    ddx = candidate_pt - pores[n, :-1]

                    # Use shortest image distance if that point lies across periodic boundary
                    ddxx = np.logical_not(periodic) * ddx + periodic * np.minimum(np.absolute(ddx + dx),
                                                                                  np.minimum(np.absolute(ddx),
                                                                                             np.absolute(ddx - dx)))
                    # print(candidate_pt, pores[n, :-1], ddxx, np.sqrt(np.sum(ddxx**2)), np.sqrt(np.sum(ddxx**2)) < spacing)

                    if np.sum(ddxx**2) < spacing**2:
                        accept = False
                        break  # Once a conflicting point is found no need to search the remaining cells

            if accept:
                # Add pore to list of pores and seed points plus mark grid cell as occupied
                grid = grid.place(tuple(candidate_cell), pores.shape[0])
                new_pore = np.append(candidate_pt, spacing/2)
                pores = np.append(pores, np.reshape(new_pore, (1, dim+1)), axis=0)
                seed_pts = np.append(seed_pts, np.reshape(candidate_pt, (1, dim)), axis=0)
            else:
                num_rejected = num_rejected + 1

    return pores


# Creates pores which some minimum spacing via random sequential adsorption (RSA)
# N = number of desired particles
# spacing = minimum separation between spheres
# r = function alias for sampling radii (To be implemented, right now only monodisperse)
# seed = random seed for RSA
# x0 = corner with minimum x, y, and z coordinates
# dx = width, height, and length of domain
# trials = number of rejected trials before terminating RSA
def rsa(N, spacing, **kwargs):
    # Setup parameters
    seed = kwargs.get('seed', np.random.randint(0, 1e8))
    dim = kwargs.get('dim', 3)
    x0 = np.array(kwargs.get('x0', np.zeros(dim)))
    dx = np.array(kwargs.get('dx', np.ones(dim)))
    trials = kwargs.get('trials', 30)
    periodic = kwargs.get('periodic', [True for k in range(dim)])

    # Grid cell size if determined given longest cell size diagonal is radius of spacing
    grid = Grid(spacing / math.sqrt(dim), dim=dim)

    pores = np.empty((0, dim + 1))

    num_rejected = 0
    np.random.seed(seed)
    while pores.shape[0] < N and num_rejected <= trials:
        x = np.empty(dim + 1)
        x[dim] = spacing/2
        for k in range(dim):
            if periodic[k]:
                x[k] = np.random.uniform(0.0, dx[k])
            else:
                x[k] = np.random.uniform(0.0+spacing/2, dx[k]-spacing/2)

        candidate_cell = grid.pos2cell(x[0:-1])

        accept = True
        if grid.occupied(candidate_cell):
            accept = False
        else:
            neighbor_pores = grid.cells_in_range(candidate_cell, 2*np.ones(dim).astype(int))
            for p in neighbor_pores:
                d = x[0:dim] - pores[p, 0:-1]
                dd = d * np.logical_not(periodic) + periodic * np.minimum(np.abs(d), np.abs(d + dx))
                dd = d * np.logical_not(periodic) + periodic * np.minimum(np.abs(dd), np.abs(d - dx))

                if np.sum(np.square(dd)) < spacing ** 2:
                    accept = False
                    break

        if accept:
            grid = grid.place(candidate_cell, pores.shape[0])
            pores = np.append(pores, np.reshape(x, (1, dim + 1)), axis=0)
            num_rejected = 0
        else:
            num_rejected = num_rejected + 1

    for i in range(len(pores)):
        pores[i, 0:dim] = pores[i, 0:dim] + x0

    return pores


#############################################
# def voronoi_volumes(points, v0):
#   # computes volumes for a voronoi tesselation of given points
#   # or returns v0 for edge points without enough data to compute.
#     v = Voronoi(points)
#     vol = np.zeros(v.npoints)
#     for i, reg_num in enumerate(v.point_region):
#         indices = v.regions[reg_num]
#         if -1 in indices: # some regions can be opened
#             vol[i] = v0
#         else:
#             vol[i] = ConvexHull(v.vertices[indices]).volume
#     return vol


# ===========================================
# END UTILITY FUNCTIONS
# ===========================================

# ===========================================
# GEOMETRY OBJECTS
# ===========================================

#############################################
class Geometry:
  """
  Base class of PFW geometry objects
  """
  @abstractmethod
  def __init__(self,
               name,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               damage=_defaultDamage,
               porosity=_defaultPorosity,
               temperature=_defaultTemperature,
               surfaceTraction=_defaultSurfaceTraction):
    self.name = name
    self.v = v if callable(v) else np.array(v) # Only need to check if callable if converting it array
    self.mat = mat
    self.group = group
    self.damage = damage
    self.particleType = particleType
    self.porosity = porosity
    self.temperature = temperature
    self.surfaceTraction = surfaceTraction if callable(surfaceTraction) else np.array(surfaceTraction)
  
  @abstractmethod
  def isInterior(self, pt):
    pass

  @abstractmethod
  def getSurfaceNormal(self, pt):
    return _defaultSurfaceNormal

  @abstractmethod
  def getSurfacePosition(self, pt):
    return _defaultSurfacePosition

  @abstractmethod
  def xMin(self):
    return -np.Inf

  @abstractmethod
  def xMax(self):
    return np.Inf


#############################################
class shell(Geometry):
  """
  Geometry object for creating a spherical shell defined by center and inner and outer radii
  """
  def __init__(self,
               name,
               x0,
               ri,
               ro,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.x0 = np.array(x0)
    self.ri = ri
    self.ro = ro

  def isInterior(self, pt, skinDepth):
    x = np.array(pt) - self.x0
    #  Is point interior?
    if ( self.ri * self.ri <= np.dot(x,x) < self.ro * self.ro ):
      #  Is point a surface?
      if ( ( np.dot(x,x) > (self.ro - skinDepth) * (self.ro - skinDepth) ) or ( np.dot(x,x) < (self.ri + skinDepth) * (self.ri + skinDepth) ) ):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1 

  def getSurfaceNormal(self,pt):
    x = np.array(pt) - self.x0
    xn = np.linalg.norm(x)

    n = x / xn 
    if np:
      n = -n
    return n

  def getSurfacePosition(self,pt):
    x = np.array(pt) - self.x0
    xn = np.linalg.norm(x)
    n = x / xn 
    di = xn - self.ri
    do = self.ro - xn

    if di < do:
      return di*-n
    else:
      return do*n

  def xMin(self):
    return self.x0[0] - self.ro

  def xMax(self):
    return self.x0[0] + self.ro


#############################################
class sphere(Geometry):
  """
  Geometry object for creating a sphere defined by center and radius
  """
  def __init__(self,
               name,
               x0,
               r,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.x0 = np.array(x0)
    self.r = r
    self.rSqr = r*r

  def isInterior(self,pt,skinDepth):
    x = np.array(pt) - self.x0
    if np.dot(x,x) < self.rSqr:
      if np.dot(x,x) > (self.r - skinDepth)*(self.r - skinDepth):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1    
  
  def getSurfaceNormal(self,pt):
    x = np.array(pt) - self.x0
    return x / np.linalg.norm( x )
  
  def getSurfacePosition(self,pt):
    x = np.array(pt) - self.x0
    mag = np.linalg.norm(x)
    return (self.r/mag-1)*x

  def xMin(self):
    return self.x0[0] - self.r

  def xMax(self):
    return self.x0[0] + self.r

#############################################
class sphericalInclusion(Geometry):
  """
  Geometry object for creating a sphere defined by center and radius, but where surface is always false (useful for inclusions)
  """
  def __init__(self,
               name,
               x0,
               r,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.x0 = np.array(x0)
    self.r = r
    self.rSqr = r*r

  def isInterior(self,pt):
    x = np.array(pt) - self.x0
    if npt.dot(x, x) < self.rSqr:
      return 0
    
    return -1

  def getSurfaceNormal(self,pt):
    return np.array([0.0,0.0,0.0])

  def getSurfacePosition(self,pt):
    return np.array([0.0,0.0,0.0])

  def xMin(self):
    return self.x0[0] - self.r

  def xMax(self):
    return self.x0[0] + self.r


#############################################
class voronoiSphericalBondedGrainComposite(Geometry):
  """
  Geometry object for creating a spherical bonded-grain composite with voronoi crystal binder inside defined by center, radius, grain size, and binder thickness
  
  # TODO fix surface flagging and implement surface position and normals
  """
  def __init__(self,
               name,
               center,
               radius,
               grainDiameter,
               binderThickness,
               shellThickness,
               porosity,
               radialBias,
               seed,
               v=_defaultVelocity,
               grainMat=0,
               binderMat=1,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               dim=3):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)    
    self.dim = dim
    self.center = np.array(center)  # center
    self.center = self.center[0:self.dim] # For 2D problem only need x and y
    self.radius = radius  # BondedGrainComposite radius inluding shell
    self.radiusSquared = radius * radius
    self.x0 = self.center - self.radius
    self.x1 = self.center + self.radius
    self.grainMat = grainMat
    self.binderMat = binderMat
    self.grainDiameter = grainDiameter
    self.binderThickness = binderThickness
    self.shellThickness = shellThickness
    self.porosity = porosity
    self.radialBias = radialBias
    self.seed = seed

    # Create evenly distributed densely packed pts to generate voronoi cells that represent grains
    self.vpts = poisson(self.grainDiameter, x0=self.x0, dx=self.x1-self.x0, seed=self.seed, dim=self.dim)
    self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
    self.npts = self.vpts.shape[0]
    self.kdt = KDTree(self.vpts, leaf_size=np.ceil(len(self.vpts) / 2), metric='euclidean')
    self.voronoi = Voronoi(self.vpts)

    # Associate a group and phase (e.g. porosity (0) or crystal (1))
    self.phase = []
    self.matDirs = []
    for i in range(self.npts):
        r = np.linalg.norm(self.vpts[i, :] - self.center)
        # Discourage porosity near the BondedGrainComposite surface to be consistent with CT imagery
        probability = 0.5 * porosity * (radialBias + 1.0) * (radialBias + 2.0) * max(0.0, (1. - r / radius) ** radialBias)
        p = (0 if (np.random.uniform(0.0, 1.0) < probability and r < radius - grainDiameter) else 1)
        self.phase.append(p)

        d = random_direction(dim=self.dim)
        if self.dim == 2:
          d = np.append(d, 0.0)

        if abs(np.dot(d, np.array([0.0,0.0,1.0]))-1) < 1e-12:
          m2 = np.cross(np.array([0.0,1.0,0.0]),d)
        else:
          m2 = np.cross(np.array([0.0,0.0,1.0]),d)
        m2 = m2 / np.linalg.norm(m2)
        m3 = np.cross(d,m2)

        self.matDirs.append(np.vstack((d, m2, m3)))

  def isInterior(self,pt, skinDepth):
    x = np.array(pt[0:self.dim])
    xx = self.center-x
    rp = np.dot(xx,xx) # np.linalg.norm(self.center-x)
    if rp < self.radiusSquared:
      # Is material in binder shell:
      if rp > (self.radius - self.shellThickness)*(self.radius - self.shellThickness):
        return True
      else:
        # material is pore IF (it is not in the binder gap AND the closest point is pore ) OR ( IF it is in the binder gap AND the closest two points are pore )

        # use the voroni cell to compute material.  We look at the distance to the two nearest points.  If the two distances are within
        # the threshold set by binderThickness, the point is on the edge of the voronoi cell, so we set it to binder otherwise it is
        # grain.
        dist, index = self.kdt.query(x.reshape(1,-1), k=2)
        d = np.abs( dist[0,1] - dist[0,0] )
        if ( d < self.binderThickness): 
          #grain boundary
          if ( self.phase[index[0,0]] == 0 and self.phase[index[0,1]] == 0):
            #boundary between pores so make it a pore
            return False
          else:
            #boundary of a grain to make it not a pore
            return True
        else:
          # not grain boundary
          if (self.phase[index[0,0]] == 1):
            #closest point is not pore, so this will be grain material
            return True
          else:
            # closest point is pore and not in boundary, so set to pore
            return False
    
    # point is outside radius of the BondedGrainComposite
    return False

  def isSurface(self,pt,skinDepth):  
    # is the particle within distance l of the surface? Includes pore surfaces and surface of neighboring grains
    # assumes the point is interior
    x = np.array(pt[0:self.dim])

    # If within skinDepth of BondedGrainComposite surface must be a surface
    xx = x - self.center
    if np.dot(xx,xx) > (self.radius - skinDepth)*(self.radius - skinDepth):
        return _defaultSurfaceFlag

    # Otherwise check if near surface of grain inside BondedGrainComposite
    # Don't need to check if voronoi cell is porosity because that should have been caught by isInterior
    dist, index = self.kdt.query(x.reshape(1, -1), k=1)
    index = index[0,0]

    # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
    minSurfaceDist = np.Inf
    ridgePts = self.voronoi.ridge_points
    for i in range(len(ridgePts)):
      p1 = ridgePts[i][0]
      p2 = ridgePts[i][1]
      if index == p1 or index == p2:
          if index == p1:
              p = p2
          else:
              p = p1

          # If neighbor cell is another grain, skip it
          if self.phase[p] == 1 and self.phase[index] == 1:
            continue

          n = self.vpts[p, :] - self.vpts[index, :]
          n = n / 2
          d = np.linalg.norm(n)
          n = n / d

          dv = x - self.vpts[index, :]
          dvc = np.dot(n, dv)  # component of points along voronoi face normal
          if dvc > 0.0:
            if self.phase[index] == 1:
              # Inside a crystal
              
              if dvc > d - skinDepth + self.binderThickness/2:
                return _defaultSurfaceFlag
            else:
              # Inside a pore
              dd = self.radius - np.linalg.norm(xx)
              if dd < dvc:
                if dd > self.radius - self.shellThickness and dd < self.radius - self.shellThickness + skinDepth:
                  return _defaultSurfaceFlag
              else:
                if dvc < d + skinDepth - self.binderThickness/2:
                  return _defaultSurfaceFlag

    return 0

  def getSurfaceNormal(self,pt):
    return np.array([0.0,0.0,0.0])

  def getSurfacePosition(self,pt):
    return np.array([0.0,0.0,0.0])

  def xMin(self):
    return self.x0[0]

  def xMax(self):
    return self.x1[0]

  def getMat(self,pt):
    # assumes the point within the object
    x = np.array(pt[0:self.dim])

    # if point is on surface, it is binder
    r = np.linalg.norm(x-self.center)
    if ( r > self.radius-self.shellThickness ):
      mat = self.binderMat
    else: 
      # use the voroni cell to compute material.  We look at the distance to the two nearest points.  If the two distances are within
      # the threshold set by binderThickness, the point is on the edge of the voronoi cell, so we set it to binder otherwise it is
      # grain.
      dist, index = self.kdt.query(x.reshape(1,-1), k=2)
      d = np.abs( dist[0,1] - dist[0,0] )
      if ( d < self.binderThickness):
        mat = self.binderMat
      else:
        mat = self.grainMat

    return mat # will error if map doesn't contain index as key

  def getMatDir(self, pt):
    x = np.array(pt[0:self.dim])
    dist, index = self.kdt.query(x.reshape(1, -1), k=1)
    index = index[0,0]

    return self.matDirs[index]


#############################################
class czSphericalBondedGrainComposite(Geometry):
    """
    Geometry object for creating a spherical BondedGrainComposite with voronoi crystals bound by cohesive zones and defined by center, radius, and grain size

    Applies to 2D and 3D BondedGrainComposites with cohesive zones defined between neighboring grains

    # TODO should take a flag to assign as single group or group offset to avoid including grains with previously defined geometry objects
    """
    def __init__(self,
                 name,
                 center,
                 radius,
                 grainDiameter,
                 porosity,
                 radialBias,
                 seed,
                 v=_defaultVelocity,
                 mat=_defaultMat,
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim=3):
        Geometry.__init__(self,
                          name,
                          v = v,
                          mat = mat,
                          group = group,
                          particleType = particleType)
        self.dim = dim
        self.center = np.array(center)  # center
        self.center = self.center[0:self.dim] # For 2D problem only need x and y
        self.radius = radius  # BondedGrainComposite radius inluding shell
        self.radiusSquared = radius * radius
        self.x0 = self.center - self.radius
        self.x1 = self.center + self.radius
        self.grainDiameter = grainDiameter
        self.porosity = porosity
        self.radialBias = radialBias
        self.seed = seed

        # Create evenly distributed densely packed pts to generate voronoi cells that represent grains
        self.vpts = poisson(self.grainDiameter, x0=self.x0, dx=self.x1-self.x0, seed=self.seed, dim=self.dim)
        self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
        self.npts = self.vpts.shape[0]
        self.kdt = KDTree(self.vpts, leaf_size=np.ceil(len(self.vpts) / 2), metric='euclidean')
        self.voronoi = Voronoi(self.vpts)

        # Associate a group and phase (e.g. porosity (0) or crystal (1))
        self.phase = []
        # self.group = []
        self.matDirs = []
        for i in range(self.npts):
            r = np.linalg.norm(self.vpts[i, :] - self.center)
            # Discourage porosity near the BondedGrainComposite surface to be consistent with CT imagery
            probability = 0.5 * porosity * (radialBias + 1.0) * (radialBias + 2.0) * max(0.0, (1. - r / radius) ** radialBias)
            p = (0 if (np.random.uniform(0.0, 1.0) < probability and r < radius - grainDiameter) else 1)
            self.phase.append(p)
            # self.group.append(i)

            d = random_direction(dim=self.dim)
            if self.dim == 2:
              d = np.append(d, 0.0)

            if abs(np.dot(d, np.array([0.0,0.0,1.0]))-1) < 1e-12:
              m2 = np.cross(np.array([0.0,1.0,0.0]),d)
            else:
              m2 = np.cross(np.array([0.0,0.0,1.0]),d)
            m2 = m2 / np.linalg.norm(m2)
            m3 = np.cross(d,m2)

            self.matDirs.append(np.vstack((d, m2, m3)))

    def isInterior(self, pt, skinDepth):
        # is the point within the object
        x = np.array(pt[0:self.dim])
        xx = x - self.center

        # Check if point is inside sphere
        if np.dot(xx, xx) < self.radiusSquared:
            # Find voronoi cell closest to point
            dist, index = self.kdt.query(x.reshape(1, -1), k=1)
            index = index[0,0]
            # If voronoi cell is not porosity
            if self.phase[index] != 0:
              surfaceFlag = 0 # Particle is interior unless otherwise determined

              # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
              minSurfaceDist = np.Inf
              ridgePts = self.voronoi.ridge_points
              for i in range(len(ridgePts)):
                  p1 = ridgePts[i][0]
                  p2 = ridgePts[i][1]
                  if index == p1 or index == p2:
                      if index == p1:
                          p = p2
                      else:
                          p = p1

                      n = self.vpts[p, :] - self.vpts[index, :]
                      n = n / 2
                      d = np.linalg.norm(n)
                      n = n / d

                      dv = x - self.vpts[index, :]
                      dvc = np.dot(n, dv)  # component of points along voronoi face normal
                      if dvc > 0.0 and dvc > d - skinDepth:
                          surfaceFlag = 3
                          break

              # If within skinDepth of BondedGrainComposite surface must be a surface
              xx = x - self.center
              if np.dot(xx,xx) > (self.radius - skinDepth)*(self.radius - skinDepth):
                if surfaceFlag != 3:
                  surfaceFlag = _defaultSurfaceFlag

              return surfaceFlag
        
        return -1
       

    def getSurfaceNormal(self, pt):
        # assumes the point is interior and a surface
        x = np.array(pt[0:self.dim])

        dist, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.Inf
        surfaceNormal = np.Inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfaceNormal = n

        # Compute distance to surface of BondedGrainComposite and check if closer
        xx = x - self.center
        xx_mag = np.linalg.norm(xx)
        surfaceDistance = self.radius - xx_mag
        if minSurfaceDist > surfaceDistance:
            surfaceNormal = xx / xx_mag

        if self.dim == 2:
          surfaceNormal = np.append(surfaceNormal, 0.0)

        return surfaceNormal

    def getSurfacePosition(self, pt):
        # assumes the point is interior and a surface
        x = np.array(pt[0:self.dim])

        dist, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.Inf
        surfacePosition = np.Inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfacePosition = surfaceDistance * n

        # Compute distance to surface of BondedGrainComposite and check if closer
        xx = x - self.center
        xx_mag = np.linalg.norm(xx)
        surfaceDistance = self.radius - xx_mag
        if minSurfaceDist > surfaceDistance:
            surfacePosition = surfaceDistance * xx / xx_mag

        if self.dim == 2:
          surfacePosition = np.append(surfacePosition, 0.0)

        return surfacePosition

    def getGroup(self, pt):
        # x = np.array(pt[0:self.dim])

        # dist, index = self.kdt.query(x.reshape(1, -1), k=1)
        # index = index[0,0]
        return  self.group # self.group[index]

    def getMatDir(self, pt):
        x = np.array(pt[0:self.dim])

        dist, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        return self.matDirs[index]

    def xMin(self):
        return self.x0[0]

    def xMax(self):
        return self.x1[0]


#############################################
class czCylindricalBondedGrainComposite(Geometry):
    """
    Geometry object for creating a cylindrical BondedGrainComposite with voronoi crystals bound by cohesive zones and defined by center, radius, and grain size

    Applies to 2D and 3D BondedGrainComposites with cohesive zones defined between neighboring grains

    # TODO should take a flag to assign as single group or group offset to avoid including grains with previously defined geometry objects
    """
    def __init__(self,
                 name,
                 x1,
                 x2,
                 radius,
                 grainDiameter,
                 porosity,
                 radialBias,
                 seed,
                 v=_defaultVelocity,
                 mat=_defaultMat,
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim=3):
        Geometry.__init__(self,
                          name,
                          v = v,
                          mat = mat,
                          group = group,
                          particleType = particleType)
        self.dim = dim
        self.x1 = np.array(x0)
        self.x2 = np.array(x1)
        self.a = self.x1-self.x0
        self.a = self.a/np.linalg.norm(self.a)
        self.center = 0.5*(self.x0+self.x1)  # center
        self.radius = radius  # BondedGrainComposite radius inluding shell
        self.radiusSquared = radius * radius

        self.grainDiameter = grainDiameter
        self.porosity = porosity
        self.radialBias = radialBias
        self.seed = seed

        # Create evenly distributed densely packed pts to generate voronoi cells that represent grains
        self.vpts = poisson(self.grainDiameter, x0=self.x0, dx=self.x1-self.x0, seed=self.seed, dim=self.dim)
        self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
        self.npts = self.vpts.shape[0]
        self.kdt = KDTree(self.vpts, leaf_size=np.ceil(len(self.vpts) / 2), metric='euclidean')
        self.voronoi = Voronoi(self.vpts)

        # Associate a group and phase (e.g. porosity (0) or crystal (1))
        self.phase = []
        self.matDirs = []
        for i in range(self.npts):
            r = np.linalg.norm(self.vpts[i, :] - self.center)
            # Discourage porosity near the BondedGrainComposite surface to be consistent with CT imagery
            probability = 0.5 * porosity * (radialBias + 1.0) * (radialBias + 2.0) * max(0.0, (1. - r / radius) ** radialBias)
            p = (0 if (np.random.uniform(0.0, 1.0) < probability and r < radius - grainDiameter) else 1)
            self.phase.append(p)

            d = random_direction(dim=self.dim)
            if self.dim == 2:
              d = np.append(d, 0.0)

            if abs(np.dot(d, np.array([0.0,0.0,1.0]))-1) < 1e-12:
              m2 = np.cross(np.array([0.0,1.0,0.0]),d)
            else:
              m2 = np.cross(np.array([0.0,0.0,1.0]),d)
            m2 = m2 / np.linalg.norm(m2)
            m3 = np.cross(d,m2)

            self.matDirs.append(np.vstack((d, m2, m3)))

    def insideCylinder(self, pt, l):
      # THIS doesn't seem quite right, FIXME
      # Seems like it's a check of surface flag not a check of interiorness?
      x = np.array(pt) - self.x1
      z = np.dot(x, self.axis)                                   # z-coordinate of test point
      r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point

      if (z<l or z>(self.h-l)) or (r > self.r-l):
        return True
      
      return False

    def isInterior(self, pt, skinDepth):
      # is the point within the object
      x = np.array(pt[0:self.dim])

      # Check if point is inside sphere
      if self.insideCylinder(pt, 0.0):
          # Find voronoi cell closest to point
          dist, index = self.kdt.query(x.reshape(1, -1), k=1)
          index = index[0,0]

          # If voronoi cell is not porosity
          if self.phase[index] != 0:
            surfaceFlag = 0 # Particle is interior unless otherwise determined

            # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
            minSurfaceDist = np.Inf
            ridgePts = self.voronoi.ridge_points
            for i in range(len(ridgePts)):
                p1 = ridgePts[i][0]
                p2 = ridgePts[i][1]
                if index == p1 or index == p2:
                    if index == p1:
                        p = p2
                    else:
                        p = p1

                    n = self.vpts[p, :] - self.vpts[index, :]
                    n = n / 2
                    d = np.linalg.norm(n)
                    n = n / d

                    dv = x - self.vpts[index, :]
                    dvc = np.dot(n, dv)  # component of points along voronoi face normal
                    if dvc > 0.0 and dvc > d - skinDepth:
                        surfaceFlag = 3
                        break

            # If within skinDepth of BondedGrainComposite surface must be a surface
            if isInterior(pt, skinDepth):
              if surfaceFlag != 3:
                surfaceFlag = _defaultSurfaceFlag

            return surfaceFlag
      
      return -1      

    def getSurfaceNormal(self, pt):
        # assumes the point is interior and a surface
        x = np.array(pt[0:self.dim])

        dist, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.Inf
        surfaceNormal = np.Inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfaceNormal = n

        # Compute distance to surface of BondedGrainComposite and check if closer
        xx = x - self.center
        xx_mag = np.linalg.norm(xx)
        surfaceDistance = self.radius - xx_mag
        if minSurfaceDist > surfaceDistance:
            surfaceNormal = xx / xx_mag

        if self.dim == 2:
          surfaceNormal = np.append(surfaceNormal, 0.0)

        return surfaceNormal

    def getSurfacePosition(self, pt):
      # assumes the point is interior and a surface
      x = np.array(pt[0:self.dim])

      dist, index = self.kdt.query(x.reshape(1, -1), k=1)
      index = index[0,0]

      minSurfaceDist = np.Inf
      surfacePosition = np.Inf*np.ones((1,self.dim))
      ridgePts = self.voronoi.ridge_points
      for i in range(len(ridgePts)):
          p1 = ridgePts[i][0]
          p2 = ridgePts[i][1]
          if index == p1 or index == p2:
              if index == p1:
                  p = p2
              else:
                  p = p1
              n = self.vpts[p, :] - self.vpts[index, :]
              n = n / 2
              d = np.linalg.norm(n)
              n = n / d

              dv = (x - self.vpts[index, :])
              surfaceDistance = (d - np.dot(dv, n))
              if minSurfaceDist > surfaceDistance:
                  minSurfaceDist = surfaceDistance
                  surfacePosition = surfaceDistance * n

      x = np.array(pt)-self.x1
      z = np.dot(x,self.axis)  # z-coordinate of test point
      xr = x - z*self.axis
      r = np.linalg.norm( xr )  # r coordinate of test point

      dist_from_wall = self.r-r
      dist_from_top = self.h-z
      dist_from_bot = z

      min_surf_dist = min(dist_from_wall, min(dist_from_top, dist_from_bot))

      if min_surf_dist == dist_from_top:
        return dist_from_top * self.axis
      
      if min_surf_dist == dist_from_bot:
        return dist_from_bot * -self.axis

      if min_surf_dist == dist_from_wall:
        return (self.r/r - 1) * xr

      # Compute distance to surface of BondedGrainComposite and check if closer
      xx = x - self.center
      xx_mag = np.linalg.norm(xx)
      surfaceDistance = self.radius - xx_mag
      if minSurfaceDist > surfaceDistance:
          surfacePosition = surfaceDistance * xx / xx_mag

      if self.dim == 2:
        surfacePosition = np.append(surfacePosition, 0.0)

      return surfacePosition

    def getGroup(self, pt):
        return  self.group

    def getMatDir(self, pt):
        x = np.array(pt[0:self.dim])

        dist, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        return self.matDirs[index]


#############################################
class czBondedGrainComposite(Geometry):
    """
    Geometry object for creating a BondedGrainComposite (box) with voronoi crystals bound by cohesive zones and defined by minimum corner, maximum corner, and grain size

    Applies to 2D and 3D BondedGrainComposites with cohesive zones defined between neighboring grains

    # TODO should take a flag to assign as single group or group offset to avoid including grains with previously defined geometry objects
    """
    def __init__(self,
                 name,
                 xmin,
                 xmax,
                 grainDiameter,
                 porosity,
                 radialBias,
                 seed,
                 v=_defaultVelocity,
                 mat=_defaultMat,
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 bondedSurfaceFraction=1.0,
                 dim=3):
        Geometry.__init__(self,
                          name,
                          v = v,
                          mat = mat,
                          group = group,
                          particleType = particleType)
        self.dim = dim
        self.xmin = np.array(xmin[0:self.dim])
        self.xmax = np.array(xmax[0:self.dim])
        self.dx = self.xmax - self.xmin
        self.center = 0.5 * ( self.xmin + self.xmax )
        self.grainDiameter = grainDiameter
        self.porosity = porosity
        self.radialBias = radialBias
        self.seed = seed

        self.bondedSurfaceFraction = bondedSurfaceFraction

        # Create evenly distributed densely packed pts to generate voronoi cells that represent grains
        self.vpts = poisson(self.grainDiameter, x0=self.xmin, dx=self.dx, seed=self.seed, dim=self.dim)
        self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
        self.npts = self.vpts.shape[0]
        self.kdt = KDTree(self.vpts, leaf_size=np.ceil(len(self.vpts) / 2), metric='euclidean')
        self.voronoi = Voronoi(self.vpts)

        # Associate a group and phase (e.g. porosity (0) or crystal (1))
        self.phase = []
        self.matDirs = []
        radius = np.min((self.xmax-self.xmin)*0.5)
        for i in range(self.npts):
            r = np.linalg.norm(self.vpts[i, :] - self.center)

            # Discourage porosity near the BondedGrainComposite surface to be consistent with CT imagery
            probability = 0.5 * porosity * (radialBias + 1.0) * (radialBias + 2.0) * max(0.0, (1. - r / radius) ** radialBias)
            p = (0 if (np.random.uniform(0.0, 1.0) < probability and r < radius - grainDiameter) else 1)
            self.phase.append(p)

            d = random_direction(dim=self.dim)
            if self.dim == 2:
              d = np.append(d, 0.0)

            if abs(np.dot(d, np.array([0.0,0.0,1.0]))-1) < 1e-12:
              m2 = np.cross(np.array([0.0,1.0,0.0]),d)
            else:
              m2 = np.cross(np.array([0.0,0.0,1.0]),d)
            m2 = m2 / np.linalg.norm(m2)
            m3 = np.cross(d,m2)

            self.matDirs.append(np.vstack((d, m2, m3)))

        # This should define a surface flag for each cell face based on the specified fraction of bonded surfaces.
        self.ridgePtFlags = [ (3 if ( np.random.random() <= self.bondedSurfaceFraction ) else 2 ) for p in self.voronoi.ridge_points ]

    def isInterior(self, pt, skinDepth):
        # is the point within the object
        x = np.array(pt[0:self.dim])

        # Check if point is inside bounding
        if inside_box(x, self.xmin, self.dx, [False for d in range(self.dim)]):
          # Find voronoi cell closest to point
          dist, index = self.kdt.query(x.reshape(1, -1), k=1)
          index = index[0,0]
          
          # If voronoi cell is not porosity
          if self.phase[index] != 0:
            # Iterate over all ridge points and check if it is skinDepth distance from voronoi face
            minSurfaceDist = np.Inf
            ridgePts = self.voronoi.ridge_points
            for i in range(len(ridgePts)):
                p1 = ridgePts[i][0]
                p2 = ridgePts[i][1]
                if index == p1 or index == p2:
                    if index == p1:
                        p = p2
                    else:
                        p = p1

                    n = self.vpts[p, :] - self.vpts[index, :]
                    n = n / 2
                    d = np.linalg.norm(n)
                    n = n / d

                    dv = x - self.vpts[index, :]
                    dvc = np.dot(n, dv)  # component of points along voronoi face normal
                    if dvc > 0.0 and dvc > d - skinDepth:
                        return self.ridgePtFlags[i]
          return 0

        return -1
       

    def getSurfaceNormal(self, pt):
        # assumes the point is interior and a surface
        x = np.array(pt[0:self.dim])

        dist, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.Inf
        surfaceNormal = np.Inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfaceNormal = n

        if self.dim == 2:
          surfaceNormal = np.append(surfaceNormal, np.array([0.0]))

        return surfaceNormal

    def getSurfacePosition(self, pt):
        # assumes the point is interior and a surface
        x = np.array(pt[0:self.dim])

        dist, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        minSurfaceDist = np.Inf
        surfacePosition = np.Inf*np.ones((1,self.dim))
        ridgePts = self.voronoi.ridge_points
        for i in range(len(ridgePts)):
            p1 = ridgePts[i][0]
            p2 = ridgePts[i][1]
            if index == p1 or index == p2:
                if index == p1:
                    p = p2
                else:
                    p = p1
                n = self.vpts[p, :] - self.vpts[index, :]
                n = n / 2
                d = np.linalg.norm(n)
                n = n / d

                dv = (x - self.vpts[index, :])
                surfaceDistance = (d - np.dot(dv, n))
                if minSurfaceDist > surfaceDistance:
                    minSurfaceDist = surfaceDistance
                    surfacePosition = surfaceDistance * n

        if self.dim == 2:
          surfacePosition = np.append(surfacePosition, np.array([0.0]))

        return surfacePosition

    def getGroup(self, pt):
        return  self.group

    def getMatDir(self, pt):
        x = np.array(pt[0:self.dim])

        dist, index = self.kdt.query(x.reshape(1, -1), k=1)
        index = index[0,0]

        return self.matDirs[index]

    def xMin(self):
        return self.xmin[0]

    def xMax(self):
        return self.xmax[0]


#############################################
class ellipsoid(Geometry):
  """
  Geometry object for creating a grid aligned ellipsoid defined by center and three lengths
  """
  def __init__(self,
               name,
               x0,
               a,
               b,
               c,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.x0 = np.array(x0)
    self.a = a
    self.b = b
    self.c = c
    self.tolerance = 1e-6

  def isInterior(self,pt,skinDepth):
    x = pt[0] - self.x0[0]
    y = pt[1] - self.x0[1]
    z = pt[2] - self.x0[2]

    if (x/self.a)*(x/self.a) + (y/self.b)*(y/self.b) + (z/self.c)*(z/self.c) < 1:
      a = max( skinDepth/100, self.a - skinDepth)
      b = max( skinDepth/100, self.b - skinDepth)
      c = max( skinDepth/100, self.c - skinDepth)

      if ( (x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c) > 1 ):
        return _defaultSurfaceFlag
      
      return 0
    
    return -1

  def getSurfaceNormal(self,pt):
    p = self.getSurfacePosition(pt)
    n = np.divide(p, np.power(np.array([self.a,self.b,self.c]),2))
    return n / np.linalg.norm(n)

  def getSurfacePosition(self,pt):
    x = np.array(pt) - self.x0
    L = max(self.a,max(self.b,self.c))
    p1 = np.array([0.0,0.0,0.0])
    p2 = x
    n = p2 / np.linalg.norm(p2)
    p3 = p2 + L*n

    # Bisection method (numerical method to find surface point on ellipsoid)
    max_iterations = 100
    for i in range(max_iterations):
      g2 = np.sum(np.divide(np.power(p2,2),np.power(np.array([self.a,self.b,self.c]),2)))
      if abs(g2 - 1.0) < self.tolerance:
        return p2 - x

      if g2 < 1:
        p1 = p2
      else:
        p3 = p2

      p2 = (p1+p3)/2
    
    return None # Should through error    

  def xMin(self):
    return -self.a

  def xMax(self):
    return self.a
  

#############################################
class crystal(Geometry):
  """
  Generates faceted crystal analog in 3D defined by a height and min and max face offsets
  """
  def __init__(self,
               name,
               center,
               axis,
               height,
               rmin,
               rmax,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.center = center
    self.axis = axis
    self.height = height
    self.rmin = rmin
    self.rmax = rmax

    n0 = np.array(self.axis)
    n0 = n0 / np.linalg.norm(n0)

    n1 = np.array([3/(5*np.sqrt(2)), 2*np.sqrt(2)/5, 1/np.sqrt(2)]) 
    if np.all(n0 == n1):
      n1 = np.array([3/(5*np.sqrt(2)), 1/np.sqrt(2), 2*np.sqrt(2)/5])

    n1 = np.dot(np.identity(3)-np.tensordot(n0, n0, axes=0), n1)
    n1 = n1 / np.linalg.norm(n1)
    n2 = np.cross(n0, n1)
    self.matDir = np.vstack((n0, n1, n2))

    faceNormals = []

    # Face set 1
    for theta in np.linspace(0.0, 5*math.pi/3, 6): #math.pi/3):
      nn = np.cos(theta)*n1 + np.sin(theta)*n2
      nn = nn / np.linalg.norm(nn)
      faceNormals.append(nn)

    # Face set 2
    for theta in np.linspace(math.pi/6, 3*math.pi/2, 3): # 2*math.pi/3):
      nn = n0 + np.cos(theta)*n1 + np.sin(theta)*n2
      nn = nn / np.linalg.norm(nn)
      faceNormals.append(nn)

    # Face set 3
    for theta in np.linspace(math.pi/6, 3*math.pi/2, 3): #2*math.pi/3):
      nn = -n0 + np.cos(theta)*n1 + np.sin(theta)*n2
      nn = nn / np.linalg.norm(nn)
      faceNormals.append(nn)

    faceNormals.append(n0)
    faceNormals.append(-n0)

    self.faceNormals = faceNormals

    self.ci = np.vstack((np.zeros((6,3)), 
                         np.tile(self.height/2*n0,(3,1)),
                         np.tile(-self.height/2*n0,(3,1)),
                         np.zeros((2,3))))
    self.ri = (self.rmax-self.rmin)*np.random.rand(14,1)+self.rmin
    self.ri[-2:] = self.ri[-2:] + self.height

  def isInterior(self, pt, skinDepth):
    pc = np.array(pt) - self.center
    surfaceFlag = 0 # Surface flag defaults to interior unless otherwise determined
    for i in range(len(self.faceNormals)):
      dist2Surf = np.dot(self.faceNormals[i], pc-self.ci[i,:])
      if dist2Surf > self.ri[i]:
        return -1
      
      if dist2Surf > self.ri[i] - skinDepth:
        surfaceFlag = _defaultSurfaceFlag

    return surfaceFlag

  def getSurfaceNormal(self,pt):
    pc = np.array(pt) - self.center
    min_distance_2_surf = np.Inf
    surface_normal = np.empty((1,3))
    for i in range(len(self.faceNormals)):
      m = self.ri[i] - np.dot(self.faceNormals[i], pc-self.ci[i,:])
      if m < min_distance_2_surf:
        min_distance_2_surf = m
        surface_normal = self.faceNormals[i]

    return surface_normal

  def getSurfacePosition(self,pt):
    pc = np.array(pt) - self.center
    min_distance_2_surf = np.Inf
    surface_normal = np.empty((1,3))
    for i in range(len(self.faceNormals)):
      m = self.ri[i] - np.dot(self.faceNormals[i], pc-self.ci[i,:])
      if m < min_distance_2_surf:
        min_distance_2_surf = m
        surface_normal = self.faceNormals[i]

    return min_distance_2_surf * surface_normal

  def xMin(self):
    return -rmax

  def xMax(self):   
    return rmax


#############################################
class indentor(Geometry):
  """
  Geometry object for creating a grid-aligned indentor defined by angle and number of facets'

  note: alpha=65.3 deg for berkovich with nFaces=3
  note: alpha=68 deg for vickers with nFaces=4
  """
  def __init__(self,
               name,
               x0,
               nFaces,
               alpha,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.x0 = np.array(x0)         # tip coordinate, assumes identor is aligned with z-axis with tip in -z direction
    self.nFaces = nFaces # number of facets
    self.alpha = np.radians(alpha)   # facet angle, input in degrees, stored in radians, 
    
    # create set of unit vectors orthogonal to each facet.
    theta = 2.*np.pi/nFaces # angle between normals about indentor axis
    self.n = []
    k=0.
    for i in range(0,nFaces):
      k=k+1.
      ny = -np.sin(self.alpha)
      r = np.cos(self.alpha)
      nx = r*np.cos(theta*(k-0.5))
      nz = r*np.sin(theta*(k-0.5))

      # nz = -np.sin(self.alpha)
      # r = np.cos(self.alpha)
      # nx = r*np.cos(theta*(k-0.5))
      # ny = r*np.sin(theta*(k-0.5))

      (self.n).append(np.array([nx,ny,nz]))

  def isInterior(self, pt, skinDepth):  
    x = np.array(pt)-self.x0
    surfaceFlag = 0
    for i in range(0,self.nFaces):
      dist2Surf = x.dot(self.n[i]) 
      if dist2Surf > 0.0:
        return -1

      if dist2Surf > -skinDepth:
        surfaceFlag = _defaultSurfaceFlag
    
    return surfaceFlag

  def getSurfaceNormal(self,pt):
    x = np.array(pt)-self.x0
    min_distance_2_surf = np.Inf
    surface_normal = np.empty((1,3))
    for i in range(0,self.nFaces):
      m =  -x.dot(self.n[i])
      if m < min_distance_2_surf:
        min_distance_2_surf = m
        surface_normal = self.n[i]

    return surface_normal

  def getSurfacePosition(self,pt):
    x = np.array(pt)-self.x0
    min_distance_2_surf = np.Inf
    surface_normal = np.empty((1,3))
    for i in range(0,self.nFaces):
      m =  -x.dot(self.n[i])
      if m < min_distance_2_surf:
        min_distance_2_surf = m
        surface_normal = self.n[i]

    return min_distance_2_surf * surface_normal


#############################################
class spherical_indenter(Geometry):
  """
  Geometry object for creating a grid-aligned spherical tip indentor defined by radius and angle
  """
  def __init__(self,
               name,
               x0,
               radius,
               alpha,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               dim=3):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.x0 = np.array(x0)          # tip coordinate (bottom of spherical tip), assumes identor is aligned with y-axis with tip in -y direction
    self.radius = radius
    self.rSqr = radius * radius
    self.alpha = np.radians(alpha)   # tip half angle, input in degrees, stored in radians
    self.dim = dim

    self.g = self.x0[1] + self.radius * (1.0 - np.sin(self.alpha) )
    self.vtip = self.g - np.cos(self.alpha) * self.radius / np.tan(self.alpha) # height of virtual tip of cone shape
    
  def isInterior(self, pt, skinDepth):  
    x = np.array(pt)-np.array([self.x0[0], 0.0, self.x0[2]])
    if self.dim == 2:
      x[2] = 0.0

    surfaceFlag = 0
    if x[1] >= self.g:
      x[1] = x[1] - self.vtip
      rr = x[0] * x[0] + x[2] * x[2]
      if rr < (x[1] * np.tan(self.alpha))**2:
        if rr >= (x[1] * np.tan(self.alpha) - skinDepth * np.cos(self.alpha))**2:
          return _defaultSurfaceFlag
        else:
          return 0
    else:
      x[1] = x[1] - self.radius - self.x0[1]
      if np.dot(x,x) < self.rSqr:
        if np.dot(x,x) > (self.radius - skinDepth) * (self.radius - skinDepth):
          return _defaultSurfaceFlag
        else:
          return 0
    
    return -1    

  def getSurfaceNormal(self,pt):
    x = np.array(pt)-np.array([self.x0[0], 0.0, self.x0[2]])
    if self.dim == 2:
      x[2] = 0.0

    if x[1] >= self.g:
      m = np.sqrt( x[0] * x[0] + x[2] * x[2]  )
      h = m * np.tan(self.alpha)

      n = np.array([x[0], -h, x[2]])
      surfaceNormal = n / np.linalg.norm(n)
    else:
      x[1] = x[1] - self.radius - self.x0[1]
      surfaceNormal = x / np.linalg.norm(x)
    
    return surfaceNormal

  def getSurfacePosition(self,pt):
    x = np.array(pt)-np.array([self.x0[0], 0.0, self.x0[2]])
    if self.dim == 2:
      x[2] = 0.0

    if x[1] >= self.g:
      m = np.sqrt( x[0] * x[0] + x[2] * x[2]  )
      h = m * np.tan(self.alpha)

      n = np.array([x[0], -h, x[2]])
      n = n / np.linalg.norm(n)
      x[1] =  x[1] - self.vtip
      surfacePosition = -np.dot(n, x) * n
    else:
      x[1] = x[1] - ( self.radius + self.x0[1] )
      xmag = np.linalg.norm(x)
      n = x / xmag
      surfacePosition = (self.radius - xmag) * n 
    
    return surfacePosition


#############################################
class box(Geometry):
    """
    Geometry object for creating a grid-aligned box defined by two corners.
    """
    def __init__(self, 
                 name, 
                 x0, 
                 x1, 
                 v=_defaultVelocity, 
                 mat=_defaultMat, 
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim=3, 
                 flaggedSurfaces=[True, True, True, True, True, True]):
      Geometry.__init__(self,
                        name,
                        v = v,
                        mat = mat,
                        group = group,
                        particleType = particleType)
      self.dim = dim
      self.x0 = np.array(x0[0:self.dim])
      self.x1 = np.array(x1[0:self.dim])
      self.xmin = np.minimum(self.x0, self.x1)
      self.xmax = np.maximum(self.x0, self.x1)
      self.flaggedSurfaces=np.array(flaggedSurfaces[0:2*self.dim])

    def isInterior(self, pt, skinDepth):
      x = np.array(pt[0:self.dim])
      if np.all( np.logical_and( x >= self.xmin, x < self.xmax) ):
        s = np.hstack((x <= self.xmin + skinDepth, x >= self.xmax - skinDepth))
        if  np.any( np.logical_and( s, self.flaggedSurfaces ) ):
          return _defaultSurfaceFlag
        else:
          return 0
      
      return -1

    def getSurfaceNormal(self, pt):
      x = np.array(pt[0:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.Inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      minI = np.argmin(np.absolute(dx))

      d = minI % self.dim
      s = -1 if int(math.floor(minI / self.dim) == 0) else 1

      surfaceNormal = np.array([0.0, 0.0, 0.0])
      surfaceNormal[d] = s
      return surfaceNormal

    def getSurfacePosition(self, pt):
      x = np.array(pt[0:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.Inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      minI = np.argmin(np.absolute(dx))

      d = minI % self.dim
      s = -1 if int(math.floor(minI / self.dim) == 0) else 1

      surfacePosition = np.array([0.0, 0.0, 0.0])
      surfacePosition[d] = dx[minI]
      return surfacePosition

      # x = np.array(pt[0:self.dim])
      # dx = np.vstack((self.x0 - x, self.x1 - x))

      # dxI = np.argmin(np.absolute(dx), axis=0)
      # dxMin = []
      # for d in range(self.dim):
      #   dxMin.append(dx[dxI[d]][d])
      # dxMin = np.array(dxMin)

      # minI = np.argmin(np.absolute(dxMin))

      # surfacePos = np.array([0.0, 0.0, 0.0])
      # surfacePos[minI] = dx[dxI[minI]][minI]

      # return surfacePos

    def xMin(self):
      return self.xmin[0]

    def xMax(self):
      return self.xmax[0]

#############################################
class notchedBar(Geometry):
    """
    Grid-aligned box defined by two corners with an edge notch in the +y face, having a 45degree
    angle and specified depth

    NOTE, TODO: Surface normals around notch have not been defined.
    """
    def __init__(self, 
                 name, 
                 x0, # corner 1 
                 x1, # corner 2
                 h, # notch depth
                 v=_defaultVelocity, 
                 mat=_defaultMat, 
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim=3, 
                 flaggedSurfaces=[True, True, True, True, True, True]):
      Geometry.__init__(self,
                        name,
                        v = v,
                        mat = mat,
                        group = group,
                        particleType = particleType)
      self.dim = dim
      self.x0 = np.array(x0[0:self.dim])
      self.x1 = np.array(x1[0:self.dim])
      self.h = h
      self.xmin = np.minimum(self.x0, self.x1)
      self.xmax = np.maximum(self.x0, self.x1)
      self.flaggedSurfaces=np.array(flaggedSurfaces[0:2*self.dim])

    def isInterior(self, pt, skinDepth):
      x = np.array(pt[0:self.dim])
      if ( np.all( np.logical_and( x >= self.xmin, x < self.xmax) ) and ( x[1] < self.xmax[1] - self.h + np.abs(x[0] - 0.5*( self.xmin[0] + self.xmax[0] ) ) ) ):    
        s = np.hstack((x <= self.xmin + skinDepth, x >= self.xmax - skinDepth))
        if  np.any( np.logical_and( s, self.flaggedSurfaces ) ):
          return _defaultSurfaceFlag
        else:
          return 0
      
      return -1

    def getSurfaceNormal(self, pt):
      x = np.array(pt[0:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.Inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      minI = np.argmin(np.absolute(dx))

      d = minI % self.dim
      s = -1 if int(math.floor(minI / self.dim) == 0) else 1

      surfaceNormal = np.array([0.0, 0.0, 0.0])
      surfaceNormal[d] = s
      return surfaceNormal

    def getSurfacePosition(self, pt):
      x = np.array(pt[0:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.Inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      minI = np.argmin(np.absolute(dx))

      d = minI % self.dim
      s = -1 if int(math.floor(minI / self.dim) == 0) else 1

      surfacePosition = np.array([0.0, 0.0, 0.0])
      surfacePosition[d] = dx[minI]
      return surfacePosition

      # x = np.array(pt[0:self.dim])
      # dx = np.vstack((self.x0 - x, self.x1 - x))

      # dxI = np.argmin(np.absolute(dx), axis=0)
      # dxMin = []
      # for d in range(self.dim):
      #   dxMin.append(dx[dxI[d]][d])
      # dxMin = np.array(dxMin)

      # minI = np.argmin(np.absolute(dxMin))

      # surfacePos = np.array([0.0, 0.0, 0.0])
      # surfacePos[minI] = dx[dxI[minI]][minI]

      # return surfacePos

    def xMin(self):
      return self.xmin[0]

    def xMax(self):
      return self.xmax[0]

#############################################
class box2(Geometry):
    """
    Geometry object for creating a grid-aligned box defined by two corners, mixed normals at corners for contact testing
    """
    def __init__(self, 
                 name, 
                 x0, 
                 x1, 
                 v=_defaultVelocity, 
                 mat=_defaultMat, 
                 group=_defaultGroup,
                 particleType=_defaultParticleType,
                 dim=3, 
                 flaggedSurfaces=[True, True, True, True, True, True]):
      Geometry.__init__(self,
                        name,
                        v = v,
                        mat = mat,
                        group = group,
                        particleType = particleType)
      self.dim = dim
      self.x0 = np.array(x0[0:self.dim])
      self.x1 = np.array(x1[0:self.dim])
      self.xmin = np.minimum(self.x0, self.x1)
      self.xmax = np.maximum(self.x0, self.x1)
      self.flaggedSurfaces=np.array(flaggedSurfaces[0:2*self.dim])

    def isInterior(self, pt, skinDepth):
      x = np.array(pt[0:self.dim])
      if np.all( np.logical_and( x >= self.xmin, x < self.xmax) ):
        s = np.hstack((x <= self.xmin + skinDepth, x >= self.xmax - skinDepth))
        if  np.any( np.logical_and( s, self.flaggedSurfaces ) ):
          return _defaultSurfaceFlag
        else:
          return 0
      
      return -1

    def getSurfaceNormal(self, pt):
      x = np.array(pt[0:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.Inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      dx_min = np.min(np.absolute(dx))
      nearest = ( np.absolute( dx ) - dx_min )< 1e-6
      

      surfaceNormal = np.array( [dx[0] if nearest[0] else 0.0 + dx[2] if nearest[2] else 0.0, dx[1] if nearest[1] else 0.0 + dx[3] if nearest[3] else 0.0, 0.0])
      # print(dx, dx_min, nearest, surfaceNormal, np.linalg.norm(surfaceNormal))
      # surfaceNormal = np.array( [nearest[0] * -dx[0] + nearest[3] * dx[3], nearest[1] * -dx[1] + nearest[4] * dx[4], nearest[2] * -dx[2] + nearest[5] * dx[5]])
      # if self.dim == 2:
      #   surfaceNormal[2] = 0.0
      return surfaceNormal / np.linalg.norm(surfaceNormal)

    def getSurfacePosition(self, pt):
      x = np.array(pt[0:self.dim])
      m = np.zeros((2*self.dim))
      m[np.logical_not(self.flaggedSurfaces)] = np.Inf
      dx = np.hstack((self.xmin - x, self.xmax - x)) + m
      dx_min = np.min(np.absolute(dx))
      nearest = ( np.absolute( dx ) - dx_min )< 1e-2

      surfaceNormal = np.array( [dx[0] if nearest[0] else 0.0 + dx[2] if nearest[2] else 0.0, dx[1] if nearest[1] else 0.0 + dx[3] if nearest[3] else 0.0, 0.0])
      # if self.dim == 2:
      #   surfaceNormal[2] = 0.0
      return surfaceNormal

    def xMin(self):
      return self.xmin[0]

    def xMax(self):
      return self.xmax[0]


#############################################
class polygon:
  """
  Geometry object for creating a polygon described by ordered vertices.'
  """
  def __init__(self,
               name,
               plist,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType,
               flaggedSurfaces=None):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.plist = np.array(plist)
    self.flaggedSurfaces = flaggedSurfaces
    if self.flaggedSurfaces is None:
      self.flaggedSurfaces = [True for i in range(self.plist.shape[0])]

  def ccw(self,x1,x2,x3):
    return (x3[1]-x1[1])*(x2[0]-x1[0]) - (x2[1]-x1[1])*(x3[0]-x1[0]) > -10**-16

  def intersect(self,A,B,C,D):
    return self.ccw(A,C,D) != self.ccw(B,C,D) and self.ccw(A,B,C) != self.ccw(A,B,D)

  def isInside(self,vertices,point):
    v_arr = np.array(vertices)
    xmin = min(v_arr[:,0])
    ymin = min(v_arr[:,1])
    xmax = max(v_arr[:,0])
    ymax = max(v_arr[:,1])
    dx = xmax - xmin
    dy = ymax - ymin
    outside = [xmin,ymin]
    outside[0] -= dx
    outside[1] -= dy
    nv = v_arr.shape[0]
    nIntersections = 0
    for i in range(nv):
      p1 = vertices[i]
      p2 = vertices[np.mod(i+1,nv)]
      if(self.intersect(p1,p2,outside,point)):
        nIntersections += 1
    if np.mod(nIntersections,2)==0:
      return False
    else:
      return True

  def isInterior(self, point, skinDepth):
    x = np.array(point)
    if(self.isInside(self.plist,x)):
      vertices = self.plist
      nv = vertices.shape[0]

      nearestI = -1
      nearestEdgeD = np.Inf
      for i in range(nv):
        if not self.flaggedSurfaces[i]:
          continue
        
        j = np.mod(i+1,nv)
        v = x - vertices[i,:]
        w = vertices[j,:] - vertices[i,:]

        wNorm = np.linalg.norm(w)
        w = w / wNorm
        dw = np.dot(v,w)
        if dw >= 0.0 and dw <= wNorm:
          v = v - dw * w
          d = np.linalg.norm(v)
          if d < nearestEdgeD:
            nearestEdgeD = d
            nearestI = i

      if nearestEdgeD < skinDepth:
        return _defaultSurfaceFlag
      
      return 0
    
    return -1

  def getSurfaceNormal(self,pt):
    x = np.array(pt)
    # Assumes point is already internal
    # Find the nearest edge and use it's normal
    vertices = self.plist
    nv = vertices.shape[0]

    nearestEdgeD = np.Inf
    surfaceNormal = np.empty((1,3))
    for i in range(nv):
      j = np.mod(i+1,nv)
      v = x - vertices[i,:]
      w = vertices[j,:] - vertices[i,:]

      wNorm = np.linalg.norm(w)
      w = w / wNorm
      dw = np.dot(v,w)
      if dw >= 0.0 and dw <= wNorm:
        v = v - dw * w
        d = np.linalg.norm(v)
        if d < nearestEdgeD:
          nearestEdgeD = d
          surfaceNormal = -v / d

    return surfaceNormal
  
  def getSurfacePosition(self,pt):
    x = np.array(pt)
    # Assumes point is already internal
    # Find the nearest edge and vector to surface
    vertices = self.plist
    nv = vertices.shape[0]
    
    nearestEdgeD = np.Inf
    surfaceNormal = np.empty((1,3))
    for i in range(nv):
      j = np.mod(i+1,nv)
      v = x - vertices[i,:]
      w = vertices[j,:] - vertices[i,:]

      wNorm = np.linalg.norm(w)
      w = w / wNorm
      dw = np.dot(v,w)
      if dw >= 0.0 and dw <= wNorm:
        v = v - dw * w
        d = np.linalg.norm(v)
        if d < nearestEdgeD:
          nearestEdgeD = d
          surfacePosition = -v

    return surfacePosition

  def xMin(self):
    arr = np.array(self.plist)
    xmin = min(arr[:,0])
    return xmin

  def xMax(self):
    arr = np.array(self.plist)
    xmax = max(arr[:,0])
    return xmax


#############################################
class foam(Geometry):
  """
  Geometry object for creating a grid-aligned box defined by two corners with spherical pores.
  
  This fills a box with spherical pores defined by the array pores
  assumes this is a numpy array [[r,x,y,z],[r,x,y,z],...]
  The surface flags should be correctly set and the searching of pores
  for isInterior uses a KD-tree. which should be faster than
  using particle file writer with each pore as its own object.
  """
  def __init__(self,
               name,
               x0,
               x1,
               pores,
               v=_defaultVelocity, 
               mat=_defaultMat, 
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)
    self.pores = pores

    pts = []
    for p in pores:
      pts.append([ p[1], p[2], p[3] ])
    
    # neighbor list for points:
    self.kdt = KDTree(pts, leaf_size=np.ceil(len(pts)/2), metric='euclidean')

  def isInterior(self,pt,skinDepth):   
    x = np.array(pt)
    
    if np.all(np.logical_and(x >= self.x0, x < self.x1)):
      dist, index = self.kdt.query(x.reshape(1,-1), k=5)
      surfaceFlag = 0 # Interior unless otherwise determined
      for i in index[0]:
        p=self.pores[i]
        rSqr = ( x[0] - p[1] )**2 + ( x[1] - p[2] )**2 + ( x[2] - p[3] )**2
        if rSqr < p[0]**2:
          return -1 # Inside a pore
        
        if rSqr < ( p[0] + skinDepth )**2:
          surfaceFlag = _defaultSurfaceFlag # near internal pore surface
          
      return surfaceFlag
    
    return -1

  def xMin(self):
    return self.x0[0]

  def xMax(self):
    return self.x1[0]


#############################################
class twoFieldBox(Geometry):
  """
  Geometry object for creating a grid-aligned box defined by two corners with random assignment of group 1 or group 2.
  """
  def __init__(self,
               name,
               x0,
               x1,
               v=_defaultVelocity,
               mat=_defaultMat,
               group1=0,
               group2=1,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)

  def isInterior(self, pt, skinDepth):
    x = np.array(pt)
    # is the point within the object?
    if np.all( np.logical_and(x > self.x0, x < self.x1) ):
      # is point on the surface?
      if np.any( np.logical_or( x <= self.x0 + skinDepth, x >= self.x1- skinDepth) ):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1

  def xMin(self):
    return min(self.x0[0], self.X1[0])

  def xMax(self):
    return max(self.x0[0], self.X1[0])


#############################################
class cylinder(Geometry):
  """
  Geometry object for creating a cylinder defined by two points and radius
  """
  def __init__(self,
               name,
               x1,
               x2,
               r,
               ri=0.,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.x1 = np.array(x1)
    self.x2 = np.array(x2)
    self.r = r
    self.ri = ri #inner radius, defaults to 0 for solid cylinder
    self.center = self.x2-self.x1 # Temporary to get surface normals of 2D disks

    self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    self.axis = (self.x2-self.x1)/self.h    

    # for now we will default to cylinder having axis-aligned fiber direction,
    # may make more general later, but this allows for backwards compatibility

    if abs(np.dot(self.axis, np.array([0.0,0.0,1.0]))-1) < 1e-12:
      m2 = np.cross(np.array([0.0,1.0,0.0]),self.axis)
    else:
      m2 = np.cross(np.array([0.0,0.0,1.0]),self.axis)
    m2 = m2 / np.linalg.norm(m2)
    m3 = np.cross(self.axis,m2) # TODO double check this

    self.matDir = np.vstack((self.axis, m2, m3))

  def isInterior(self, pt, skinDepth):    
    x = np.array(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point

    if (z >= 0 and z < self.h) and (r < self.r) and (self.ri==0 or r > self.ri):
      if ( (z<skinDepth) or (z>self.h-skinDepth) or ( r > self.r-skinDepth) or (self.ri==0 or r < self.ri+skinDepth) ):
        return _defaultSurfaceFlag
      else:
        return 0    
    return -1

  def getSurfaceNormal(self,pt):
    x = np.array(pt)-self.x1
    z = np.dot(x, self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point

    dist_from_inner_wall = r - self.ri
    dist_from_outer_wall = self.r - r
    dist_from_top = self.h-z
    dist_from_bot = z

    if (self.ri > 0):
        min_surf_dist = min( [ dist_from_inner_wall, dist_from_outer_wall, dist_from_top, dist_from_bot ] )
    else:
        min_surf_dist = min( [ dist_from_outer_wall, dist_from_top, dist_from_bot ] )

    if min_surf_dist == dist_from_top:
      return self.axis

    if min_surf_dist == dist_from_bot:
      return -self.axis

    if min_surf_dist == dist_from_outer_wall:
      return xr / r

    if min_surf_dist == dist_from_inner_wall:
      return -xr / r

  def getSurfacePosition(self,pt):
    x = np.array(pt)-self.x1
    z = np.dot(x,self.axis)  # z-coordinate of test point
    xr = x - z*self.axis
    r = np.linalg.norm( xr )  # r coordinate of test point

    dist_from_inner_wall = r - self.ri
    dist_from_outer_wall = self.r - r
    dist_from_top = self.h-z
    dist_from_bot = z

    if (self.ri > 0):
        min_surf_dist = min( [ dist_from_inner_wall, dist_from_outer_wall, dist_from_top, dist_from_bot ] )
    else:
        min_surf_dist = min( [ dist_from_outer_wall, dist_from_top, dist_from_bot ] )

    if min_surf_dist == dist_from_top:
      return dist_from_top * self.axis
    
    if min_surf_dist == dist_from_bot:
      return dist_from_bot * -self.axis

    if min_surf_dist == dist_from_outer_wall:
      return (self.r/r - 1) * xr

    if min_surf_dist == dist_from_inner_wall:
      return (self.ri/r - 1) * xr

  def xMin(self):
    return min(self.x1[0], self.x2[0])-self.r

  def xMax(self):
    return max(self.x1[0], self.x2[0])+self.r


# #############################################
# class expandingRing(Geometry):
#   """
#   Geometry object for creating a cylinder defined by two points and radius
#   """
#   def __init__(self,
#                name,
#                x1,
#                x2,
#                innerRadius,
#                outerRadius,
#                radialVelocity,
#                v=_defaultVelocity,
#                mat=_defaultMat,
#                group=_defaultGroup,
#                particleType=_defaultParticleType):
#     Geometry.__init__(self,
#                       name,
#                       v = v,
#                       mat = mat,
#                       group = group,
#                       particleType = particleType)
#     self.x1 = np.array(x1)
#     self.x2 = np.array(x2)
#     self.innerRadius = innerRadius
#     self.outerRadius = outerRadius

#     self.center = self.x2-self.x1 # Temporary to get surface normals of 2D disks
#     self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
#     self.axis = (self.x2-self.x1)/self.h    

#     # for now we will default to cylinder having axis-aligned fiber direction,
#     # may make more general later, but this allows for backwards compatibility

#     if abs(np.dot(self.axis, np.array([0.0,0.0,1.0]))-1) < 1e-12:
#       m2 = np.cross(np.array([0.0,1.0,0.0]),self.axis)
#     else:
#       m2 = np.cross(np.array([0.0,0.0,1.0]),self.axis)
#     m2 = m2 / np.linalg.norm(m2)
#     m3 = np.cross(self.axis,m2) # TODO double check this

#     self.matDir = np.vstack((self.axis, m2, m3))

#   def isInterior(self, pt, skinDepth):    
#     x = np.array(pt)-self.x1
#     z = np.dot(x,self.axis)  # z-coordinate of test point
#     r = np.linalg.norm( x - z*self.axis )  # r coordinate of test point

#     if (z >= 0 and z < self.h) and (r < self.outerRadius) and (r > self.innerRadius):
#       if ( (z<skinDepth) or (z>self.h-skinDepth) or ( r > self.outerRadius-skinDepth) or (r < self.innerRadius+skinDepth) ):
#         return _defaultSurfaceFlag
#       else:
#         return 0
    
#     return -1

#   def getSurfaceNormal(self,pt):
#     x = np.array(pt)-self.x1
#     z = np.dot(x, self.axis)  # z-coordinate of test point
#     xr = x - z*self.axis
#     r = np.linalg.norm( xr )  # r coordinate of test point

#     dist_from_inner_wall = r - self.innerRadius
#     dist_from_outer_wall = self.outerRadius - r
#     dist_from_top = self.h-z
#     dist_from_bot = z

#     min_surf_dist = min( [ dist_from_inner_wall, dist_from_outer_wall, dist_from_top, dist_from_bot ] )

#     if min_surf_dist == dist_from_top:
#       return self.axis

#     if min_surf_dist == dist_from_bot:
#       return -self.axis

#     if min_surf_dist == dist_from_inner_wall:
#       return -xr / r

#     if min_surf_dist == dist_from_outer_wall:
#       return xr / r

#   def getSurfacePosition(self,pt):
#     x = np.array(pt)-self.x1
#     z = np.dot(x,self.axis)  # z-coordinate of test point
#     xr = x - z*self.axis
#     r = np.linalg.norm( xr )  # r coordinate of test point

#     dist_from_inner_wall = r - self.innerRadius
#     dist_from_outer_wall = self.outerRadius - r
#     dist_from_top = self.h-z
#     dist_from_bot = z

#     min_surf_dist = min( [ dist_from_inner_wall, dist_from_outer_wall, dist_from_top, dist_from_bot ] )

#     if min_surf_dist == dist_from_top:
#       return dist_from_top * self.axis
    
#     if min_surf_dist == dist_from_bot:
#       return dist_from_bot * -self.axis

#     if min_surf_dist == dist_from_inner_wall:
#       return (self.innerRadius/r - 1) * xr

#     if min_surf_dist == dist_from_outer_wall:
#       return (self.outerRadius/r - 1) * xr

#   def getVelocity(self,pt):
#     x = np.array(pt)-self.x1
#     z = np.dot(x,self.axis)  # z-coordinate of test point
#     xr = x - z*self.axis
#     r = np.linalg.norm( xr )  # r coordinate of test point

#     if ( r > 1.e-12):
#       er = xr/r
#       v = self.v + self.radialVelocity*er
#     else:
#       v = self.v

#     return v

#   def xMin(self):
#     return min(self.x1[0], self.x2[0]) - max(self.innerRadius,self.outerRadius) 

#   def xMax(self):
#     return max(self.x1[0], self.x2[0]) + max(self.innerRadius,self.outerRadius) 

#############################################
class whiskers(Geometry):
  """
  Geometry object for creating a cylinder defined by two points and radius
  """
  def __init__(self,
               name,
               x1,
               x2,
               r,
               numWhiskers,
               volFracWhiskers,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.x1 = np.array(x1)
    self.x2 = np.array(x2)
    self.r = r
    self.nw = numWhiskers
    self.VF = volFracWhiskers

    self.h = np.linalg.norm(self.x2-self.x1) # height of cylinder axis
    self.e1 = (self.x2-self.x1)/self.h       # unit vector along cylinder axis

    # random vectors orthogonal to e1 and each other
    self.e2 = np.random.rand(3)
    self.e2 = self.e2 - np.dot(self.e2,self.e1)*self.e1
    self.e2 = self.e2/np.linalg.norm(self.e2)
    self.e3 = np.cross(self.e1,self.e2)

  def isInterior(self, pt, skinDepth):    
    x = np.array(pt)
    z = np.dot(x-self.x1,self.e1)                                   # z-coordinate of test point
    pMinusX1InPlane = (x-self.x1) - np.dot(x-self.x1,self.e1)*self.e1
    r = np.linalg.norm( pMinusX1InPlane )  # r coordinate of test point
  
    ex = np.dot( self.e2,pMinusX1InPlane )
    ey = np.dot( self.e3,pMinusX1InPlane )

    theta = np.arctan2(ex,ey)
    dtheta = 2.0*np.pi/self.nw

    if ( ( (z>=0 and z<self.h) and (r < self.r) ) and ( (theta%dtheta)/dtheta < self.VF) ):
      if ( (z<skinDepth or z>(self.h-skinDepth)) or (r > self.r-skinDepth) ):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1

  def xMin(self):
    return min(self.x1[0], self.x2[0])-self.r

  def xMax(self):
    return max(self.x1[0], self.x2[0])+self.r


#############################################
class toroid(Geometry):
  """
  Geometry object for creating a circular toroid defined by a point, direction ring radius and revolved circle radius
  """
  def __init__(self,
               name,
               x0,
               n,
               r1,
               r2,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.x0 = np.array(x0) # center point
    self.n = np.array(n) # axes unit vector
    self.n = self.n / np.linalg.norm(self.n)
    self.r1 = r1         # ring radius
    self.r2 = r2         # revolved circle radius

  def isInterior(self,pt,skinDepth):    
    x = np.array(pt) - self.x0
    z = x.dot(self.n)
    r = np.linalg.norm(x - z*self.n)

    if (np.sqrt( (r-self.r1)*(r-self.r1) + z*z ) < self.r2):
      if (np.sqrt( (r-self.r1)*(r-self.r1) + z*z ) > self.r2 - skinDepth):
        return _defaultSurfaceFlag
      else:
        return 0
    
    return -1

  def getSurfaceNormal(self,pt):
    x = np.array(pt) - self.x0
    z = x.dot(self.n)
    xr = x - z*self.n
    r = np.linalg.norm(xr)
    xrr = x - self.r1 * xr / r
    xrr_norm = np.linalg.norm(xrr)

    return xrr / xrr_norm

  def getSurfacePosition(self,pt):
    x = np.array(pt) - self.x0
    z = x.dot(self.n)
    xr = x - z*self.n
    r = np.linalg.norm(xr)
    xrr = x - self.r1 * xr / r 
    xrr_norm = np.linalg.norm( xrr )

    return ( self.r2 / xrr_norm - 1) * xrr

  def xMin(self):
    return self.x0[0]-self.r1-self.r2

  def xMax(self):
    return self.x0[0]+self.r1+self.r2


#############################################
class spinodal(Geometry):
  """
  Geometry object for generating arbitrary tailorable periodic spinodals microstructures using wavelet method

  Because N must be large it is too expensive to generate spinodal directly. Instead distance field is generated and interpolated to 
  generate spinodal microstructure

  # rho = target relative density
  # seed  = random seed for generating wave vectors and phases
  # num_waves = number of waves should be large (1000-100000)
  # a = cell size specified as a scalar or vector [ax, ay, az]
  # shell = boolean determines if geometry is a shell or solid spinodal
  """
  def __init__(self, 
               name, 
               rho, 
               seed, 
               shell=False, 
               grid_size=[100,100,100], 
               a=[0.33333,0.33333,0.33333], 
               num_waves=1000,
               v=_defaultVelocity, 
               mat=_defaultMat, 
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.rho = rho
    self.seed = seed
    self.a = np.array(a)
    self.cellSize = np.array([a,a,a])
    self.shell = shell
    self.num_waves = num_waves
    self.offset = np.array([0.0,0.0,0.0])
    if self.shell:
      rho_s = self.rho / 2 + 0.5
    else:
      rho_s = self.rho
    self.level_set = np.sqrt( 2.0 ) * scipy.special.erfinv(2.0*rho_s-1.0)
    
    print(self.rho, rho_s, self.level_set)

    self.A = np.sqrt( 2.0 ) / self.num_waves #np.sqrt( 2.0 / self.num_waves )  # Prefactor for level set function (avoid computing for every point)

    self.phase = []
    self.q = []
    for n in range(self.num_waves):
      self.q.append(np.divide(2*np.pi*random_direction(), self.a).tolist())
      self.phase.append(np.random.uniform(0.0, 2*np.pi))
    
    self.phase  = np.array(self.phase)
    self.q = np.array(self.q)
    print(self.q.shape, self.phase.shape)

  def fphi(self, u):
    c = self.A*np.sum(np.cos(np.inner(self.q, u[:,np.newaxis].T)+self.phase))
    # print(c, self.A, np.sum(np.cos(np.inner(self.q, u[:,np.newaxis].T)+self.phase)), self.level_set, np.inner(self.q, u[:,np.newaxis].T).shape)
    
    return c

  def fdphi(self, u):
    c = 0
    for n in range(self.num_waves):
      c = c - self.kappa * np.sin(self.kappa * np.dot(self.k[n], u) + self.phase[n])
    
    return (1.41421356237 * c / self.num_waves) 

  # def pt2coords(self, pt):
  #   return 2*np.pi*np.divide(np.array(pt) - self.offset, self.cellSize) 

  def vec2Surf(self, pt):
    # Currently seem to work with gyroid may need to be adjusted later
    tolerance = 0.01
    maxIter = 1000
    alpha = 0.1
    
    c = self.level_set
    p = pt
    phi = self.fphi(p)
    dphi = self.fdphi(p)
    # print("i:", 0, ", p:",p,", c:", self.level_set,", phi:", phi, ", dphi", dphi, ", de:", phi-c)
    for i in range(maxIter):
      s = 1
      if phi <= 0:
        c = -abs(c)
        s = 0
      dpt = (-1)**(1+s)*alpha*abs(phi - c)*dphi
      phi = self.fphi(p)
      dphi = self.fdphi(p)
      p = p + dpt
      # print("i:", i, ", p:",p,", c:", self.level_set,", phi:", phi, ", dphi", dphi, ", de:", phi-c)
      if abs(phi - c) < tolerance:
        # print("Converged!")
        return p - pt
      
    print('Did not converge!')

  def isInterior(self, pt, skinDepth):
    # u = self.pt2coords(pt)
    u = np.array(pt)

    # Use level set method
    if self.shell:
      if abs(self.fphi(u)) < self.level_set:
        return True
    else:
      # print(u, self.fphi(u), self.level_set)
      if self.fphi(u) < self.level_set:
        return True

    return False

  def getSurfaceNormal(self, pt):
    # u = self.pt2coords(pt)
    # v = self.vec2Surf(u)
    # s = self.fdphi(u+v)
    # if self.fphi(u) <= 0:
    #   s = -s
    # return s / np.linalg.norm(s)
    return np.array([0.0,0.0,0.0])

  def getSurfacePosition(self, pt):
    # u = self.pt2coords(pt)
    # v = self.vec2Surf(u)
    # return np.multiply(v, self.cellSize)/(2*np.pi)
    return np.array([0.0,0.0,0.0])


#############################################
# Crook implementation
class tpms(Geometry):
  """
  Geometry object for creating grid aligned triply periodic minimal surface

  # tpms_type = type of triply periodic minimal surface (TPMS) or spinodal
  # rho = target relative density
  # cellSize = array of size 3 that determines cellSize in x, y, z optionally also a single scalar
  # offset = array of size 3 to offset the cell of the TPMS
  """
  def __init__(self,
               name,
               tpms_type,
               rho=0.5,
               cellSize=[1.0,1.0,1.0],
               offset=[0.0,0.0,0.0],
               v=_defaultVelocity, 
               mat=_defaultMat, 
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.tpms_type = tpms_type
    self.rho = rho
    self.cellSize= np.array(cellSize)
    self.offset = np.array(offset)

    # When python3 is updated to 3.10, we can use match instead of ifs
    if self.tpms_type == "gyroid":
      self.level_set = np.interp(self.rho, _gyroid_rhoVsC[:,0], _gyroid_rhoVsC[:,1])
      self.phi = lambda u : np.sin(u[0])*np.cos(u[1])+np.sin(u[1])*np.cos(u[2])+np.sin(u[2])*np.cos(u[0])
      self.dphi = lambda u : np.array([np.cos(u[0])*np.cos(u[1])-np.sin(u[2])*np.sin(u[0]),
                                                -np.sin(u[0])*np.sin(u[1])+np.cos(u[1])*np.cos(u[2]),
                                                -np.sin(u[1])*np.sin(u[2])+np.cos(u[2])*np.cos(u[0])])
      return

    if self.tpms_type == "shwarz diamond":
      self.level_set = np.interp(self.rho, _schwarzDiamond_rhoVsC[:,0], _gyroid_rhoVsC[:,1])
      self.phi = lambda u : np.cos(u[0])*np.cos(u[1])*np.cos(u[2]) - np.sin(u[0])*np.sin(u[1])*np.sin(u[2])
      self.dphi = lambda u : np.array([-np.sin(u[0])*np.cos(u[1])*np.cos(u[2]) - np.cos(u[0])*np.sin(u[1])*np.sin(u[2]),
                                        -np.cos(u[0])*np.sin(u[1])*np.cos(u[2]) - np.sin(u[0])*np.cos(u[1])*np.sin(u[2]),
                                        -np.cos(u[0])*np.cos(u[1])*np.sin(u[2]) - np.sin(u[0])*np.sin(u[1])*np.cos(u[2])])
      return

    if self.tpms_type == "schwarz primitive":
      self.level_set = np.interp(self.rho, _schwarzPrimitive_rhoVsC[:,0], _gyroid_rhoVsC[:,1])
      self.phi = lambda u : np.cos(u[0]) + np.cos(u[1]) + np.cos(u[2])
      self.dphi = lambda u : np.array([-np.sin(u[0]),
                                        -np.sin(u[1]),
                                        -np.sin(u[2])])
      return
    
    if self.tpms_type == "schwarz hexagonal":
      print("schwarz hexagonal tpms not currently implemented!")
      return
    
    if self.tpms_type == "schwarz clp":
      print("schwarz hexagonal tpms not currently implemented!")
      return
        
    print('No TPMS of type specified!')

  def pt2coords(self, pt):
    return 2*np.pi*np.divide(np.array(pt) - self.offset, self.cellSize) 

  def vec2Surf(self, pt):
    # Currently seem to work with gyroid may need to be adjusted later
    tolerance = 0.01
    maxIter = 1000
    alpha = 0.1
    
    c = self.level_set
    p = pt
    phi = self.phi(p)
    dphi = self.dphi(p)
    # print("i:", 0, ", p:",p,", c:", self.level_set,", phi:", phi, ", dphi", dphi, ", de:", phi-c)
    for i in range(maxIter):
      s = 1
      if phi <= 0:
        c = -abs(c)
        s = 0
      dpt = (-1)**(1+s)*alpha*abs(phi - c)*dphi
      phi = self.phi(p)
      dphi = self.dphi(p)
      p = p + dpt
      # print("i:", i, ", p:",p,", c:", self.level_set,", phi:", phi, ", dphi", dphi, ", de:", phi-c)
      if abs(phi - c) < tolerance:
        # print("Converged!")
        return p - pt
      
    print('Did not converge!')

  def isInterior(self, pt, skinDepth):
    u = self.pt2coords(pt)

    # Use level set method
    if abs(self.phi(u)) < self.level_set:
      v = self.vec2Surf(u)
      vMag = np.linalg.norm(np.multiply(v, self.cellSize))/(2*np.pi)

      if vMag < skinDepth:
        return _defaultSurfaceFlag

      return 0

    return -1

  def getSurfaceNormal(self, pt):
    u = self.pt2coords(pt)
    v = self.vec2Surf(u)
    s = self.dphi(u+v)
    if self.phi(u) <= 0:
      s = -s
    return s / np.linalg.norm(s)

  def getSurfacePosition(self, pt):
    u = self.pt2coords(pt)
    v = self.vec2Surf(u)
    return np.multiply(v, self.cellSize)/(2*np.pi)


#############################################
class fill(Geometry):
  """
  Geometry object for creating particles to fill the background.
  """
  def __init__(self,
               name,
               mat=_defaultMat,
               group=_defaultGroup):
    print('Warning: Fill should only be used with untransformed objects.')
    self.v = np.array([0.,0.,0.])
    self.name = name
    self.mat = mat
    self.group = group

  def isInterior(self, pt, skinDepth):
    return 0

  def getSurfaceNormal(self, pt):
    return np.array([0.0,0.0,0.0])

  def getSurfacePosition(self, pt):
    return np.array([0.0,0.0,0.0])


#############################################
class VCCTL(Geometry):
  """
  Geometry object for creating a 3D object from a VCCTL voxelized dataset

  data        text file with header stripped, containing one line per voxel and integer
              values indicating phase
  ni, nj, nk  number of voxels in the x,y, and z directions.
  x0, x1      coordinates of the -,+ corners of the object in the domain
  map         dictionary of mappings from index to mat#: dict([(1, 2), (3, 4)])
  mat         default material
  group       contact group

  """
  def __init__(self,
               name,
               data,
               ni,
               nj,
               nk,
               x0,
               x1,
               map,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):    
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.data = data
    self.ni = ni
    self.nj = nj
    self.nk = nk
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0, x1)
    self.x1 = np.maximum(x0, x1)

  def isInterior(self,pt,skinDepth):
    x = np.array(pt)
    if np.all( np.logical_and( x > self.x0, x < self.x1 ) ):
      i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      k = int( np.floor( self.nk*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
      i=max(min(self.ni-1,i),0)
      j=max(min(self.nj-1,j),0)
      k=max(min(self.nk-1,k),0)

      n = i*self.nj*self.nk + j*self.nk + k

      index = int(self.data[n])
      mat = self.map.get(index)

      if (mat >-1):
        if np.any( np.logical_or( x <= self.x0 + skinDepth, x >= self.x1- skinDepth) ):
          return _defaultSurfaceFlag
        
        # check if point is next to porosity      
        imin = int( max(0,i-1) )
        imax = int( min(self.ni-1,i+1) )
        jmin = int( max(0,j-1) )
        jmax = int( min(self.nj-1,j+1) )
        kmin = int( max(0,k-1) )
        kmax = int( min(self.nk-1,k+1) )
        for i in range(imin,imax):
          for j in range(jmin,jmax):
            for k in range(kmin,kmax):
              n = i*self.nj*self.nk + j*self.nk + k
              index = self.data[n]
              mat = self.map[index] # will error if map doesn't contain index as key

              if (mat < 0):
                return _defaultSurfaceFlag

        return 0 # interior to the domain and not porosity
      
      return -1  # internal porosity
    
    return -1 # Out of domain

  def getMat(self,pt):
    # assumes the point within the object
    x = np.array(pt)

    i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
    j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
    k = int( np.floor( self.nk*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
    i=max(min(self.ni-1,i),0)
    j=max(min(self.nj-1,j),0)
    k=max(min(self.nk-1,k),0)
    n = i*self.nj*self.nk + j*self.nk + k

    index = self.data[n]
    mat = self.map[index]
    return mat # will error if map doesn't contain index as key


#############################################
class CT(Geometry):
  """
  Geometry object for creating a 3D object from a CT voxelized dataset
  
  data are generated by a mathematica script, flattening a 3D array with
  ordering [z,y,x]

  This is similar (possibly identical) to the VCCTL object, but kept separate to allow
  for specialization without compromising backwards compatibility.
  data        text file with header stripped, containing one line per voxel and integer
              values indicating phase
  ni, nj, nk  number of voxels in the x,y, and z directions.
  x0, x1      coordinates of the -,+ corners of the object in the domain
  map         dictionary of mappings from index to mat#: dict([(1, 2), (3, 4)])
  mat         default material
  group       contact group
    
  """
  def __init__(self,
               name,
               data,
               ni,
               nj,
               nk,
               x0,
               x1,
               map,
               v=_defaultVelocity,
               mat=_defaultMat,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.data = data
    self.ni = ni
    self.nj = nj
    self.nk = nk
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)

  def isInterior(self,pt, skinDepth):
    x = np.array(pt)
    
    if np.all( np.logical_and( x > self.x0, x < self.x1 ) ):
      i = int( np.floor( self.ni*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      k = int( np.floor( self.nk*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )

      i=max(min(self.ni-1,i),0)
      j=max(min(self.nj-1,j),0)
      k=max(min(self.nk-1,k),0)

      n = i*self.nj*self.nk + j*self.nk + k

      index = int(self.data[n])
      mat = self.map.get(index)

      if (mat >= 0): # interior to the domain and not porosity
        if np.any( np.logical_or( x <= self.x0 + skinDepth, x >= self.x1- skinDepth) ):
          return _defaultSurfaceFlag
    
        # check if point is next to porosity
        i = int( np.floor( self.ni*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
        j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
        k = int( np.floor( self.nk*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
        i=max(min(self.ni-1,i),0)
        j=max(min(self.nj-1,j),0)
        k=max(min(self.nk-1,k),0)
        
        imin = int( max(0,i-1) )
        imax = int( min(self.ni-1,i+1) )
        jmin = int( max(0,j-1) )
        jmax = int( min(self.nj-1,j+1) )
        kmin = int( max(0,k-1) )
        kmax = int( min(self.nk-1,k+1) )
        poreSurface = False
        for i in range(imin,imax):
          for j in range(jmin,jmax):
            for k in range(kmin,kmax):
              n = i*self.nj*self.nk + j*self.nk + k
              index = int(self.data[n])
              mat = self.map.get(index) # will error if map doesn't contain index as key
              if (mat <0):
                return _defaultSurfaceFlag
        return 0   
      
      return -1  # internal porosity
    
    return -1    # outside of domain

  def getMat(self,pt):
    # assumes the point within the object
    x = np.array(pt)

    xmin = min(self.x0[0],self.x1[0])
    xmax = max(self.x0[0],self.x1[0])
    ymin = min(self.x0[1],self.x1[1])
    ymax = max(self.x0[1],self.x1[1])
    zmin = min(self.x0[2],self.x1[2])
    zmax = max(self.x0[2],self.x1[2])

    i = int( np.floor( self.ni*(x[2]-self.x0[2])/(self.x1[2]-self.x0[2]) ) )
    j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
    k = int( np.floor( self.nk*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
    i=max(min(self.ni-1,i),0)
    j=max(min(self.nj-1,j),0)
    k=max(min(self.nk-1,k),0)
    n = i*self.nj*self.nk + j*self.nk + k

    index = int(self.data[n])
    mat = self.map.get(index)
    return mat # will error if map doesn't contain index as key


#############################################
class bitmap(Geometry):
  """
  Geometry object for creating a 2D object from an image file

  data        binary array
  ni, nj, nk  number of voxels in the x,y, and z directions.
  x0, x1      coordinates of the -,+ corners of the object in the domain
  map         dictionary of mappings from index to mat#: dict([(1, 2), (3, 4)])
  mat         default material
  group       contact group

  """
  def __init__(self,
               name,
               data,
               ni,
               nj,
               x0,
               x1,
               map,
               v=_defaultVelocity,
               group=_defaultGroup,
               particleType=_defaultParticleType):
    Geometry.__init__(self,
                      name,
                      v = v,
                      mat = mat,
                      group = group,
                      particleType = particleType)
    self.data = data
    self.ni = ni
    self.nj = nj
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)

  def isInterior(self, pt, skinDepth):
    x = np.array(pt)
    if np.all( np.logical_and(x >= self.x0 , x < self.x1) ):
      i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
      j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
      index = int(self.data[j,i])

      if (index >0):
        if np.any( np.logical_or( x <= self.x0 + skinDepth, x >= self.x1- skinDepth) ):
          return _defaultSurfaceFlag
        
        # check if point is next to porosity
        i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
        j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )

        imin = int( max(0,i-1) )
        imax = int( min(self.ni-1,i+1) )
        jmin = int( max(0,j-1) )
        jmax = int( min(self.nj-1,j+1) )

        for i in range(imin,imax):
          for j in range(jmin,jmax):
            index = int(self.data[j,i])
            if (self.data[n]==0):
              return _defaultSurfaceFlag
        
        return 0   # interior to the domain and not porosity
    
    return -1   

  def getMat(self,pt):
    x = np.array(pt)
    i = int( np.floor( self.ni*(x[0]-self.x0[0])/(self.x1[0]-self.x0[0]) ) )
    j = int( np.floor( self.nj*(x[1]-self.x0[1])/(self.x1[1]-self.x0[1]) ) )
    #n = i*self.nj + j
    index = int(self.data[j,i])

    #index = self.data[n]
    return self.map[index] # will error if map doesn't contain index as key


# ===========================================
# END GEOMETRY OBJECTS
# ===========================================


# ===========================================
# GEOMETRY WRAPPERS 
# (Modify other geometry objects as inputs)
# ===========================================

#############################################
class BaseWrapper:
  @abstractmethod
  def __init__(self,name,subObject):
    self.subObject = subObject
    self.v = subObject.v
    self.mat = subObject.mat
    self.particleType = subObject.particleType

  def isInterior(self,pt, skinDepth):
    return self.subObject.isInterior(pt, skinDepth)

  def getSurfaceNormal(self, pt):
    if hasattr(self.subObject, 'getSurfaceNormal'):
        return self.subObject.getSurfaceNormal(pt)

    if hasattr(self.subObject, 'surfaceNormal'):
        return self.subObject.surfaceNormal

    return np.array([0.0, 1.0, 0.0])
    
  def getSurfacePosition(self,pt):
    if hasattr(self.subObject, 'getSurfacePosition'):
        return self.subObject.getSurfacePosition(pt)

    if hasattr(self.subObject, 'surfacePosition'):
        return self.subObject.surfacePosition
    
    return np.array([0.0, 0.0, 0.0])

  def getVelocity(self,pt):
    if hasattr(self.subObject, 'getVelocity'):
        velocity = self.subObject.getVelocity( pt )
    else:
        velocity = self.subObject.v( self.subObject, pt ) if callable( self.subObject.v ) else self.subObject.v

    return velocity

  def getGroup(self, pt):
    if hasattr(self.subObject, 'getGroup'):
        return self.subObject.getGroup(pt)

    if hasattr(self.subObject, 'group'):
        return self.subObject.group

    return 0

  def getMatDir(self, pt):
    if hasattr(self.subObject, 'getMatDir'):
        return self.subObject.getMatDir(pt)

    if hasattr(self.subObject, 'matDir'):
        return self.subObject.matDir

    return np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])

  def getDamage(self, pt):
    if hasattr(self.subObject, 'getDamage'):
        return self.subObject.getDamage(pt)

    if hasattr(self.subObject, 'damage'):
        return self.subObject.damage

    return _defaultDamage

  def getPorosity(self, pt):
    if hasattr(self.subObject, 'getPorosity'):
        return self.subObject.getPorosity(pt)

    if hasattr(self.subObject, 'porosity'):
        return self.subObject.porosity

    return _defaultPorosity

  def getTemperature(self, pt):
    if hasattr(self.subObject, 'getTemperature'):
        return self.subObject.getTemperature(pt)

    if hasattr(self.subObject, 'temperature'):
        return self.subObject.temperature

    return _defaultPorosity

  def getSurfaceTraction(self, pt):
    if hasattr(self.subObject, 'getSurfaceTraction'):
        return self.subObject.getSurfaceTraction(pt)

    if hasattr(self.subObject, 'surfaceTraction'):
        return self.subObject.surfaceTraction

    return _defaultSurfaceTraction

  def xMin(self):
    return self.subObject.xMin()

  def xMax(self):
    return self.subObject.xMax()
  

#############################################
class materialDirectionWrapper(BaseWrapper):
  def __init__(self,name,subObject,matDir):
    BaseWrapper.__init__(self, name, subObject)

    if np.shape(matDir) == (3,):
        self.matDir = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        self.matDir[0][:] = matDir
    elif np.shape(matDir) == (3,3):
        self.matDir = np.array(matDir)
        log2file("Unsupported material direction size in this branch...")
    else:
        self.matDir = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        log2file("Unsupported material direction size...")

  def getMatDir(self, pt):
    return self.matDir


############################################
class surfaceFlagWrapper(BaseWrapper):
  def __init__(self,name,subObject,flagType):
    BaseWrapper.__init__(self, name, subObject)
    self.flagType = flagType

  def isInterior(self, pt, skinDepth):
    flag = self.subObject.isInterior(pt, skinDepth)
        
    if flag > 0:
      return self.flagType

    return flag
    

############################################
class surfaceNormalWrapper(BaseWrapper):
  def __init__(self,name,subObject,surfaceNormal):
    BaseWrapper.__init__(self, name, subObject)
    self.surfaceNormal = np.array(surfaceNormal)

  def getSurfaceNormal(self, pt):
    return self.surfaceNormal


############################################
class surfacePositionWrapper(BaseWrapper):
  def __init__(self,name,subObject,surfacePosition):
    BaseWrapper.__init__(self, name, subObject)
    self.surfacePosition = np.array(surfacePosition)

  def getSurfacePosition(self, pt):
    return self.surfacePosition


############################################
class shrinkageFlagWrapper(BaseWrapper):
  def __init__(self,name,subObject,flag):
    BaseWrapper.__init__(self, name, subObject)
    self.flag = flag


############################################
class strengthScaleWrapper(BaseWrapper):
  def __init__(self,name,subObject,strengthScale):
    BaseWrapper.__init__(self, name, subObject)
    self.strengthScale = strengthScale


############################################
class damageWrapper(BaseWrapper):
  def __init__(self,name,subObject,damage):
    BaseWrapper.__init__(self, name, subObject)
    self.damage=damage


############################################
class porosityWrapper(BaseWrapper):
  def __init__(self,name,subObject,porosity):
    BaseWrapper.__init__(self, name, subObject)
    self.porosity=porosity


############################################
class temperatureWrapper(BaseWrapper):
  def __init__(self,name,subObject,temperature):
    BaseWrapper.__init__(self, name, subObject)
    self.temperature = temperature


############################################
class voronoiWeibullBoxWrapper(BaseWrapper):
  """
  Box wrapper for another object that will be used to assign voronoi-cell weibull distribution of strength scale
  Works for 2D or 3D cases
  Box should be bigger than the subobject
  """
  def __init__(self,name,subObject,x0,x1,flawSize,weibullVolume,weibullModulus,weibullSeed,vMin,vpts=None,dim=3,randomMatDir=False):
    BaseWrapper.__init__(self, name, subObject)
    self.object = subObject
    self.dim = dim
    self.x0 = np.array(x0)
    self.x1 = np.array(x1)
    self.dx = self.x1-self.x0
    self.x0 = self.x0[0:self.dim]
    self.x1 = self.x1[0:self.dim]
    self.randomMatDir = randomMatDir
    self.seed = weibullSeed
    
    if vpts is None:
      self.vpts = poisson(flawSize, x0=self.x0, dx=self.dx[0:self.dim], seed=self.seed, dim=self.dim)
    else:
      self.vpts = vpts

    self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
    self.npts = self.vpts.shape[0]
    self.kdt = KDTree(self.vpts, leaf_size=np.ceil(len(self.vpts) / 2), metric='euclidean')
    self.voronoi = Voronoi(self.vpts)

    #Average volume (area for 2D) to assign to edge cells
    v0= np.prod(self.x1-self.x0)/self.npts
    vor = self.voronoi

    # compute volume of each voronoi cell
    vol = np.zeros(vor.npoints)
    for i, reg_num in enumerate(vor.point_region):
      indices = vor.regions[reg_num]
      if ( (-1 in indices) or ( vor.vertices[vor.regions[i]].shape[0] < 1 ) ): # some regions can be opened
        vol[i] = v0
      else:
        vol[i] = ConvexHull(vor.vertices[indices]).volume

        numInteriorVertices = 0
        numVertices = vor.vertices[vor.regions[i]].shape[0]
        
        for v in vor.vertices[indices]:
          if (dim==2):
            v = np.append(v, np.array([0.0]))
          if subObject.isInterior(v,0.0) >= 0:
            numInteriorVertices += 1

        vol[i] = vol[i]*numInteriorVertices/numVertices*(self.dx[self.dim] if self.dim==2 else 1.0)

      vol[i] = max( vol[i], vMin )

    # define the value of strength scale that will be assigned to each cell's particles.
    cellStrengthScale=[]
    for i in range(0,self.npts):
      s = ( ( weibullVolume/vol[i] )*( np.log( np.random.uniform(1e-20,1.0) )/np.log(0.5) ) )**(1.0/weibullModulus)
      cellStrengthScale.append(s)

    self.cellStrengthScale = np.array(cellStrengthScale)

    # cells can have random orientation or can inheret from subObject
    cellMatDir=[]
    if ( self.randomMatDir ): 
      for i in range(0,self.npts):
        d = random_direction()
        cellMatDir.append(d)
    self.cellMatDir = np.array(cellMatDir)

#    if ( self.randomMatDir ):  
#      def getMatDir(self,pt):
#        x = np.array(pt[0:self.dim])
#        dist, index = self.kdt.query(x.reshape(1,-1), k=1)
#        matDir = self.cellMatDir[index[0]]
#        return matDir

  def getMatDir(self,pt):
    if ( self.randomMatDir ):  
      x = np.array(pt[0:self.dim])
      dist, index = self.kdt.query(x.reshape(1,-1), k=1)
      matDir = self.cellMatDir[index[0]]
    else:
      if hasattr( self.object, 'getMatDir' ):
        matDir = object.getMatDir( pt )
      elif hasattr( self.object, 'matDir' ):
        matDir = self.object.matDir( self.object, matDir ) if callable( self.object.matDir ) else self.object.matDir
      else:
        matDir = np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]])
    return matDir

  def getStrengthScale(self,pt):
    x = np.array(pt[0:self.dim])
    dist, index = self.kdt.query(x.reshape(1,-1), k=1)
    strengthScale = self.cellStrengthScale[index[0][0]]
    return strengthScale



############################################
class pointwisePorosityWrapper(BaseWrapper):
  """
  Box wrapper for another object that will randomly define isInterior=False for a fraction of points equal
  to the specified porosity.  Not deterministic and doesn't create surface flags.  Used to introduce
  a simple explicit representation of under-resolved porosity.

  Works for 2D or 3D cases
  """
  def __init__(self,name,subObject,porosity):
    BaseWrapper.__init__(self, name, subObject)
    self.object = subObject # all properties will be inhereted from this subObject, except getDamage() 
    self.porosity = porosity

  def getMat(self,pt):
    if np.random.random() < self.porosity:
        return -1
    else:
        return self.object.getMat(pt)


############################################
class voronoiMatDirBoxWrapper(BaseWrapper):
  """
  Box wrapper for another object that will be used to assign voronoi-cell random distribution of material direction
  Works for 2D or 3D cases
  Box should be bigger than the subobject
  """
  def __init__(self,name,subObject,x0,x1,flawSize,seed=1,dim=3):
    BaseWrapper.__init__(self, name, subObject)
    self.dim = dim
    self.x0 = np.array(x0)
    self.x1 = np.array(x1)
    self.dx = self.x1-self.x0
    self.x0 = self.x0[0:self.dim]
    self.x1 = self.x1[0:self.dim]
    self.seed = seed
    
    self.vpts = poisson(flawSize, x0=self.x0, dx=self.dx[0:self.dim], seed=self.seed, dim=self.dim)
    self.vpts = self.vpts[:,0:self.dim] # Remove spacing from points
    self.npts = self.vpts.shape[0]
    self.kdt = KDTree(self.vpts, leaf_size=np.ceil(len(self.vpts) / 2), metric='euclidean')

    # define the value of strength scale that will be assigned to each cell's particles.
    cellMatDir=[]
    for i in range(0,self.npts):
      d = random_direction()
      cellMatDir.append(d)

    self.cellMatDir = np.array(cellMatDir)

  def getMatDir(self,pt):
    x = np.array(pt[0:self.dim])
    dist, index = self.kdt.query(x.reshape(1,-1), k=1)
    matDir = self.cellMatDir[index[0]]
    return matDir

#############################################
class surfaceFlagBoxWrapper(BaseWrapper):
  def __init__(self,name,x0,x1,surfaceFlag,subObject):
    BaseWrapper.__init__(self, name, subObject)
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)
    self.surfaceFlag = surfaceFlag
  
  def isInterior(self, pt, skinDepth):
    x = np.array(pt)
    flag = self.subObject.isInterior(pt, skinDepth)
    if flag > 0 and np.all( np.logical_and(x > self.x0, x < self.x1)):
      return self.surfaceFlag
    
    return flag


#############################################
class damageBoxWrapper(BaseWrapper):
  def __init__(self, name, subObject, x0, x1, damage):
    BaseWrapper.__init__(self, name, subObject)
    x0 = np.array(x0)
    x1 = np.array(x1)
    self.x0 = np.minimum(x0,x1)
    self.x1 = np.maximum(x0,x1)
    self.damage = damage

  def getDamage(self, pt):
    x = np.array(pt)
    if np.all( np.logical_and(x > self.x0, x < self.x1)):
      return self.damage
    
    return self.subObject.damage


# ===========================================
# END GEOMETRY WRAPPERS
# ===========================================

# ===========================================
# TRANSFORMS
#
# ===========================================

#############################################
# Homogeneous Transform Matrices
def translate(dx):
  return np.array([[1.0, 0.0, 0.0, dx[0]],
                   [0.0, 1.0, 0.0, dx[1]],
                   [0.0, 0.0, 1.0, dx[2]],
                   [0.0, 0.0, 0.0, 1.0]])


def scale(ds):
  if len(ds) == 1:
    ds = np.array([ds,ds,ds])
  
  return np.array([[ds[0], 0.0, 0.0, 0.0],
                   [0.0, ds[1], 0.0, 0.0],
                   [0.0, 0.0, ds[2], 0.0],
                   [0.0, 0.0,   0.0, 1.0]])


# a0 = normal of reflection plane
# x0 = center of reflection
def reflect(a0, x0=np.array([0.0,0.0,0.0])):
  x0 = np.array(x0)
  a0 = np.array(a0)
  a0 = a0 / np.linalg.norm(a0)
  reflection = np.array([[1-2*a0[0]**2,     -2*a0[0]*a0[1],  -2*a0[0]*a0[2], 0.0],
                          [-2*a0[0]*a0[1], 1-2*a0[1]**2,     -2*a0[1]*a0[2], 0.0],
                          [-2*a0[0]*a0[2],  -2*a0[1]*a0[2], 1-2*a0[2]**2,    0.0],
                          [0.0,                        0.0,          0.0,    1.0]])
  return np.matmul(translate(x0), np.matmul(reflection, translate(-x0)))


# a0 = axis of rotation
# alpha = angle of rotation (radians)
# x0 = center of rotation
def rotate(alpha,a0=np.array([0.0,0.0,1.0]),x0=np.array([0.0,0.0,0.0])):
  x0 = np.array(x0)
  a0 = np.array(a0)
  aa = np.outer(a0,a0)
  A = np.array([ [ 0,     -a0[2],      a0[1]], 
                 [ a0[2],      0,     -a0[0]], 
                 [-a0[1],  a0[0],          0] ])
  R = np.identity(4)
  R[0:3,0:3] = np.cos(alpha)*np.identity(3) + (1-np.cos(alpha))*aa + np.sin(alpha)*A
  return np.matmul(translate(x0), np.matmul(R, translate(-x0)))


#############################################
class transform(BaseWrapper):
  @abstractmethod
  def __init__(self,
               name,
               subObject,
               transform):
    BaseWrapper.__init__(self,
                         name,
                         subObject)
    self.transform = np.array(transform)
    self.inverse = np.linalg.inv(self.transform[0:3,0:3])

  def transformPoint(self, pt):
    pt = np.array(pt)
    pt = np.append(pt, 1.0)
    pt = np.matmul(self.transform, pt)
    return pt[0:3]
  
  def transformVector(self, vec):
    return np.matmul(self.inverse, np.array(vec))

  def isInterior(self, pt, skinDepth):
    return super().isInterior(self.transformPoint(pt), skinDepth)

  def getSurfaceNormal(self, pt):
    surfaceNormal = super().getSurfaceNormal(self.transformPoint(pt))
    return self.transformVector(surfaceNormal)
    
  def getSurfacePosition(self,pt):
    surfacePosition = super().getSurfacePosition(self.transformPoint(pt))
    return self.transformVector(surfacePosition)

  def getGroup(self, pt):
    return super().getGroup(self.transformPoint(pt))

  def getMatDir(self, pt):
    matDir = super().getMatDir(self.transformPoint(pt))
    return self.transformVector(matDir)

  def getDamage(self, pt):
    return super().getDamage(self.transformPoint(pt))

  def getPorosity(self, pt):
    return super().getPorosity(self.transformPoint(pt))

  def getTemperature(self, pt):
    return super().getPorosity(self.transformPoint(pt))

  def getSurfaceTraction(self, pt):
    surfaceTraction = super().getSurfaceTraction(self.transformPoint(pt))
    return self.transformVector(surfaceTraction)

  def xMin(self):
    return -np.Inf

  def xMax(self):
    return np.Inf


# #############################################
# class reflect(BaseWrapper):
#   """
#   Geometry object for creating a rotated object around point x0 and axis a0 by angle alpha
#   """
#   def __init__(self,name,subObj,x0,a0):
#     BaseWrapper.__init__(self, name, subObj)
#     self.name = name
#     self.subObj = subObj         # object to be rotated
#     self.x0 = np.array(x0)           # center of reflection
#     self.a0 = a0/np.linalg.norm(a0)  # unit vector defining reflection direction.

#   def transformPoint(self, pt):
#     pt = np.array(pt)
#     aa = np.array([
#       [self.a0[0]*self.a0[0], self.a0[0]*self.a0[1], self.a0[0]*self.a0[2]], 
#       [self.a0[1]*self.a0[0], self.a0[1]*self.a0[1], self.a0[1]*self.a0[2]], 
#       [self.a0[2]*self.a0[0], self.a0[2]*self.a0[1], self.a0[2]*self.a0[2]]
#       ])
#     return pt - 2.0*( pt - self.x0 ).dot(aa)

#   def isInterior(self,pt,skinDepth):
#     return self.subObj.isInterior(self.transformPoint(pt), skinDepth)

#   def getMat(self,pt):
#     return self.subObj.getMat(self.transformPoint(pt))


# #############################################
# class rotate(BaseWrapper):
#   """
#   Geometry object for creating a rotated object around point x0 and axis a0 by angle alpha
#   """
#   def __init__(self,name,subObject,x0,a0,alpha):
#     BaseWrapper.__init__(self, name, subObject)
#     self.x0 = np.array(x0)          # center of rotation
#     self.a0 = a0/np.linalg.norm(a0) # vector defining rotation axis
#     self.alpha = alpha              # right handed rotation (radians) about axis

#   def transformPoint(self, pt):
#     aa = np.array([
#       [self.a0[0]*self.a0[0], self.a0[0]*self.a0[1], self.a0[0]*self.a0[2]], 
#       [self.a0[1]*self.a0[0], self.a0[1]*self.a0[1], self.a0[1]*self.a0[2]], 
#       [self.a0[2]*self.a0[0], self.a0[2]*self.a0[1], self.a0[2]*self.a0[2]]
#       ])
#     A = np.array([
#       [ 0         , -self.a0[2], self.a0[1]], 
#       [ self.a0[2],           0,-self.a0[0]], 
#       [-self.a0[1],  self.a0[0],          0]
#       ])
#     R = np.cos(self.alpha)*np.identity(3) + (1-np.cos(self.alpha))*aa + np.sin(self.alpha)*A
#     return self.x0 + (R.transpose()).dot(np.array(pt)-self.x0)

#   def getMat(self,pt):
#     x = self.transformPoint(pt)
#     return self.subObject.getMat(x)

#   def isInterior(self,pt,skinDepth):
#     return self.subObject.isInterior(self.transformPoint(pt), skinDepth)

# ===========================================
# END TRANSFORMS
# ===========================================

# ===========================================
# SET OPERATIONS
# 
# ===========================================
class SetOperation(Geometry):
  @abstractmethod
  def __init__(self,name,subObjA, subObjB, defaultToA=True):
    self.subObjA = subObjA
    self.subObjB = subObjB
    self.defaultToA = defaultToA

    if self.defaultToA:
      Geometry.__init__( self,
                         name,
                         v=self.subObjA.v,
                         mat=self.subObjA.mat,
                         group=self.subObjA.group,
                         particleType=self.subObjA.particleType)
    else:
      Geometry.__init__( self,
                         name,
                         v=self.subObjB.v,
                         mat=self.subObjB.mat,
                         group=self.subObjB.group,
                         particleType=self.subObjA.particleType)

  @abstractmethod
  def isInterior(self,pt,skinDepth):
    pass

  @abstractmethod
  def getSurfaceNormal(self, pt):
    if self.defaultToA:
      return self.subObjA.getSurfaceNormal(pt)
    else:
      return self.subObjB.getSurfaceNormal(pt)

  @abstractmethod
  def getSurfacePosition(self,pt):
    if self.defaultToA:
      return self.subObjA.getSurfacePosition(pt)
    else:
      return self.subObjB.getSurfacePosition(pt)

  @abstractmethod
  def getGroup(self, pt):
    if self.defaultToA:
      return self.subObjA.getGroup(pt)
    else:
      return self.subObjB.getGroup(pt)

  @abstractmethod
  def getMatDir(self, pt):
    if self.defaultToA:
      return self.subObjA.getMatDir(pt)
    else:
      return self.subObjB.getMatDir(pt)

  # @abstractmethod
  # def getDamage(self, pt):
  #   pass

  # @abstractmethod
  # def getPorosity(self, pt):
  #   pass

  # @abstractmethod
  # def getTemperature(self, pt):
  #   pass

  # @abstractmethod
  # def getSurfaceTraction(self, pt):
  #   pass


#############################################
class union(SetOperation):
  """
  Geometry object for creating a union of two objects
  """
  def __init__(self,name,A,B):
    SetOperation.__init__(self,
                          name,
                          A,
                          B)

  def isInterior(self, pt, skinDepth):
    intA = self.subObjA.isInterior(pt, skinDepth)
    intB = self.subObjB.isInterior(pt, skinDepth)

    if ( intA >=0 or intB >=0 ):
      # Works so long as we prioritize higher surface flags e.g. cohesive > default surface > interior
      return max(surfA, surfB)

    if intB >=0:
      return surfB

    if intA >=0:
      return surfA

    return -1 # Shouldn't reach here, if it did something went wrong

  # def surfacePosition(self, pt):
  #   intA = self.subObjA.isInterior(pt)
  #   intB = self.subObjB.isInterior(pt)

  #   if ( intA and intB ):
  #     # Works so long as we prioritize higher surface flags e.g. cohesive > default surface > interior
  #     return max(surfA, surfB)

  #   if intB:
  #     return surfB

  #   if intA:
  #     return surfA

  #   return None # Shouldn't reach here, if it did something went wrong

  # def surfaceNormal(self, pt):


#############################################
class intersection(SetOperation):
  """
  Geometry object for creating an intersection of two objects
  """
  def __init__(self,name,A,B):
    SetOperation.__init__(self, name, A, B)

  def isInterior(self,pt, skinDepth):
    sA = self.subObjA.isInterior(pt,skinDepth)
    sB = self.subObjB.isInterior(pt,skinDepth)
    if sA >= 0 and sB >= 0:
      return max(sA, sB) # Should work to always return correct surface flag
    
    return -1

  def getSurfaceNormal(self, pt):
    # return _defaultSurfaceNormal
    # # Cases
    # # Both A and B define surface normals
    # # Only ones defines surface normals
    # # Both define a surface normal for the same voxel (use surface normal of )
    
    sA = self.subObjA.getSurfacePosition(pt)
    sB = self.subObjB.getSurfacePosition(pt)

    if np.dot(sA, sA) <= np.dot(sB, sB):
      return self.subObjA.getSurfaceNormal(pt)
    else:
      return self.subObjB.getSurfaceNormal(pt)

  def getSurfacePosition(self, pt):
    sA = self.subObjA.getSurfacePosition(pt)
    sB = self.subObjB.getSurfacePosition(pt)

    if np.dot(sA, sA) <= np.dot(sB, sB):
      return sA
    else:
      return sB


#############################################
class difference(SetOperation):
  """
  Geometry object for creating a difference of two objects
  """
  def __init__(self,
               name,
               A,
               B):
    SetOperation.__init__(self,
                          name,
                          A,
                          B)

  def isInterior(self,pt, skinDepth):
    if( self.A.isInterior(pt,skinDepth) ^ self.B.isInterior(pt,skinDepth) ):
      # Does not currently support surface flagging, normals or positions
      return 0
    
    return-1

  def xMin(self):
    return self.A.xMin()

  def xMax(self):
    return self.A.xMax()


# ===========================================
# END SET OPERATIONS
# ===========================================
