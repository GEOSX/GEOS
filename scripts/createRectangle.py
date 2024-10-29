import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from geosx_xml_tools import xml_formatter as xmlFormatter


def arrayToXMLArray( array ):
    stringArray = '{%s}' % ', '.join(['%1.10f' % (x) for x in array])
    return stringArray

class Rectangle:
    
    def __init__(self, *args ) -> None:
        self.name = args[0]
        if len(args) == 4:
          self.normal = np.array( compute_normal_fromAngles( args[1], args[2] ), dtype=float )
          self.origin = np.array( args[3], dtype=float )
          self.strikeAngle = args[1]
          self.dipAngle = args[2]
        elif len(args) == 3:
          self.normal = np.array( args[1], dtype=float )
          self.origin = np.array( args[2], dtype=float )
          self.strikeAngle = -1000
          self.dipAngle = -1000
        else:
            print("Error: too many arguments to constructor")
        
        self.lengthVector = np.array([0.0, 0.0, 0.0], dtype=float )
        self.widthVector = np.array([0.0, 0.0, 0.0], dtype=float )
        self.dimensions = np.array([0.0, 0.0], dtype=float )
        self.compute_tangentVectors()

    def compute_tangentVectors( self ) -> None:   
        self.lengthVector = np.array( np.cross( self.normal, np.array([0.0,1.0,0.0]) ), dtype=float )
        self.lengthVector = self.lengthVector / np.linalg.norm(self.lengthVector) 
        self.widthVector = np.array( np.cross( self.normal, self.lengthVector ), dtype=float )  
        self.widthVector = self.widthVector / np.linalg.norm(self.widthVector)

        assert np.allclose([1,1,1], [self.normal.dot(self.normal), self.lengthVector.dot(self.lengthVector), self.widthVector.dot(self.widthVector)])
        assert np.allclose([0,0,0], [self.normal.dot(self.lengthVector), self.lengthVector.dot(self.widthVector), self.normal.dot(self.widthVector)])

    def print_XML_block(self):
        
        root = ET.Element("geos")
        element = ET.SubElement(root, "Rectangle")
        element.set( "name" , self.name )
        element.set( "normal", arrayToXMLArray( self.normal ) )
        element.set( "origin", arrayToXMLArray( self.origin ) )
        element.set( "lengthVector", arrayToXMLArray( self.lengthVector ) )
        element.set( "widthVector", arrayToXMLArray( self.widthVector ) )
        element.set( "dimensions", arrayToXMLArray( self.dimensions ) )

        from xml.dom import minidom

        xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(indent="   ")
        fileName = self.name + "_" + str(self.strikeAngle)+ "_" + str(self.dipAngle) + ".xml"
        f = open( fileName, "w" )
        f.write( xmlstr )
        f.close()
        xmlFormatter.format_file( fileName )

def plotPlane(xx, yy, plane):
    # a plane is a*x+b*y+c*z+d=0
    # [a,b,c] is the normal. Thus, we have to calculate
    # d and we're set
    d = -plane.origin.dot(plane.normal)

    # calculate corresponding z
    z = (-plane.normal[0] * xx - plane.normal[1] * yy - d) * 1. /plane.normal[2]

    # plot the surface
    title = str(plane.strikeAngle) + " -  " + str(plane.dipAngle)
    fig = plt.figure( title )
    plt3d = fig.add_subplot(111, projection='3d')
    plt3d.plot_surface(xx, yy, z)
    plt3d.set_zlim3d(0, 140)

def plot(planes):

     # create x,y
    xx, yy = np.meshgrid(range(146), range(75))

    for plane in planes:
        plotPlane(xx, yy, plane)

def compute_normal_fromTriangleVertices( points ):
    V = points[1] - points[0]
    W = points[2] - points[0]
    normal = np.cross( V, W )  

    unit_normal = normal / np.linalg.norm(normal)

    print(unit_normal)
    print(np.linalg.norm(unit_normal))

    return unit_normal

def compute_normal_fromAngles( strikeAngle, dipAngle ):
    '''
    Function to compute the plane unit normal vector given the strike and dip angles
    (in degrees). The strike should lie between 0 and 360 (negative ok) and the
    dip is restricted to lie between 0 and 90 measured in the direction such that 
    when you look in the strike direction, the fault dips to your right.
    '''
    deg_to_rad = np.pi/180
    strike = strikeAngle * deg_to_rad
    dip = dipAngle * deg_to_rad
    normal = [0.0, 0.0, 0.0] 
    normal[0] = -np.sin(dip)*np.sin(strike) # x component
    normal[1] =  np.sin(dip)*np.cos(strike) # y component
    normal[2] = -np.cos(dip)                # z component

    return normal
      
def main():
    origin = [40, 37.5, 10]
    planes = []
    planes.append(Rectangle("FracturePlane", -20, 60, origin))

    # Vertex # 1; X, m:  1.2026E+03; Y, m: -8.5871E+02; Z, m:  3.3080E+02
    # Vertex # 2; X, m:  1.2038E+03; Y, m: -8.5901E+02; Z, m:  3.3003E+02
    # Vertex # 3; X, m:  1.2027E+03; Y, m: -8.5843E+02; Z, m:  3.2935E+02
    points = [ np.array([1.2026E+03, -8.5871E+02, 3.3080E+02]), 
               np.array([1.2038E+03, -8.5901E+02, 3.3003E+02]), 
               np.array([1.2027E+03, -8.5843E+02, 3.2935E+02]) ]
    
    planes.append(Rectangle("FracturePlane", compute_normal_fromTriangleVertices( points ), origin))
    plot(planes)
    for plane in planes:
        plane.print_XML_block()
    plt.show()
    
if __name__ == "__main__":
    main()

    