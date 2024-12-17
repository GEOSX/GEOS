import numpy as np
import scipy.linalg

# Analytical results for displacement jump accross the casing-cement interface
# Due to inner traction P0
def analyticDisplacementJump(P0):
	# No debonding at casing-cement interface if inner stress is compression
	ur_casing_out = 0.0

	# Casing-cement debonding if inner stress is tension
	if(P0>0):
		# Input geometric parameters
		r_casing_in = 0.1    # meter
		r_casing_out = 0.106

		# Material properties
		G_casing = 80.8e9    # Pa
		K_casing = 175e9
		lambda_casing = K_casing - 2.0 / 3.0 * G_casing

		# Rigidity of the casing annulus
		rigidity = np.array([
			[2.0 * (lambda_casing + G_casing), -2.0 * G_casing / r_casing_in / r_casing_in],
			[2.0 * (lambda_casing + G_casing), -2.0 * G_casing / r_casing_out / r_casing_out]
		])

		# Vector of force
		force = np.array([P0, 0.0])

		# Compute the coefficients describing the closed-form solutions of radial displacement
		vectorCoefficientAB = np.dot(np.linalg.inv(rigidity), force)

		coeffA_casing = vectorCoefficientAB[0]
		coeffB_casing = vectorCoefficientAB[1]

		ur_casing_out = coeffA_casing * r_casing_out + coeffB_casing / r_casing_out

	return np.abs(ur_casing_out)	# displacement jump equals to the absolute value of the radial displacement at the outer surface of the casing 
