import numpy as np

from vtkmodules.util.numpy_support import (
    vtk_to_numpy,
)

import vtk_utils

def main():
  file_name  = "/usr/workspace/cusini1/geosx/geosx_dev/GEOS_2/src/coreComponents/python/modules/geosx_mesh_doctor/test.vtu" 
  input_mesh = vtk_utils.read_mesh( file_name )
  
  nx = 300 
  ny = 300
  nz = 1

  attribute = vtk_to_numpy(input_mesh.GetCellData().GetArray("attribute"))
  attribute = attribute.reshape((ny, nx, nz))

  # Transpose the array to switch dimensions back to (nx, ny, nz)
  attribute = np.transpose(attribute, (1, 0, 2))

  for i in range(nx):
      for j in range(ny):
          for k in range(nz):
              attribute[i][j][k] = -1
              if (j > 150):
                  if(i <= 175 and i > 125 ):
                     attribute[i][j][k] = 1
              if(j <= 150 and j > 50 ):
                  if(i <= 150 and i > 125 ):
                     attribute[i][j][k] = 2

              if(j <= 150 and j > 50 ):
                  if(i <= 175 and i > 150 ):
                     attribute[i][j][k] = 3 
                                      
                       
  output_file_name  = "/usr/workspace/cusini1/geosx/geosx_dev/GEOS_2/src/coreComponents/python/modules/geosx_mesh_doctor/test_attribute.vtu" 


  vtk_utils.write_mesh(input_mesh, 
                       vtk_utils.VtkOutput(output=output_file_name, is_data_mode_binary=False))            
     

if __name__=="__main__":
    main()