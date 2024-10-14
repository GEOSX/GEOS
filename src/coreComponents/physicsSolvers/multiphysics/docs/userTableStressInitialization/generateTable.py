import math
import numpy as np


model_top = 0.0 # in m
model_base = -1000.0 # in m
num_Points = 2
z_table = np.linspace(model_base, model_top, num_Points)
pp_gradient = 0.1 # in MPa/m
sigma_xx_gradient = 0.17 # in MPa/m
sigma_yy_gradient = 0.27 # in MPa/m
sigma_zz_gradient = 0.24 # in MPa/m



PoissonRatio = 0.25
YoungModulus = 100e6 # in Pa
grainBulkModulus = 1e27 # in #Pa
bulkModulus = YoungModulus/3.0/(1.0-2.0*PoissonRatio) # in Pa
BiotCoefficient = 1.0-bulkModulus/grainBulkModulus 
porosity = 0.375 
Gravity = 9.81


pp_table = abs(z_table)*pp_gradient*100000
fluidDensity = pp_gradient*100000/Gravity
effectiveSigma_xx_table = z_table*sigma_xx_gradient*100000 + BiotCoefficient*pp_table
effectiveSigma_yy_table = z_table*sigma_yy_gradient*100000 + BiotCoefficient*pp_table
effectiveSigma_zz_table = z_table*sigma_zz_gradient*100000 + BiotCoefficient*pp_table

bulkDensity = sigma_zz_gradient*100000/Gravity
rockDensity = (bulkDensity - porosity*fluidDensity)/(1.0-porosity)
print(f'Fluid Density = {fluidDensity} kg/m3, Rock Density = {rockDensity} kg/m3')

np.savetxt('porePressure.csv', pp_table, fmt='%1.2f', delimiter=',') 

np.savetxt('effectiveSigma_xx.csv', effectiveSigma_xx_table, fmt='%1.2f', delimiter=',')

np.savetxt('effectiveSigma_yy.csv', effectiveSigma_yy_table, fmt='%1.2f', delimiter=',')  

np.savetxt('effectiveSigma_zz.csv', effectiveSigma_zz_table, fmt='%1.2f', delimiter=',')         

x_table = np.linspace(0.0, 0.0, 1)
np.savetxt('x.csv', x_table, fmt='%1.2f', delimiter=',')

y_table = np.linspace(0.0, 0.0, 1)
np.savetxt('y.csv', x_table, fmt='%1.2f', delimiter=',')

np.savetxt('z.csv', z_table, fmt='%1.2f', delimiter=',')
