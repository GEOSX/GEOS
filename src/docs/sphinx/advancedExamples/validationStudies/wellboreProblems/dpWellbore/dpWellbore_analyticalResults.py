import sys
sys.path.append('../')
import numpy as np
import elastoPlasticWellboreAnalyticalSolutions as epwAnal

#Note: this solution is developed by Chen and Abousleiman (2017) for the case with zero cohesion

def elastic_plastic_boundary_DP(sh,sv,q_ep_boundary):
	# See eq. 25 of Chen and Abousleiman 2013
	sr0 = sh - np.sqrt(sh**2 - (4.0 * sh**2 + sv**2 - 2 * sh * sv - q_ep_boundary**2) / 3.0)
	s00 = 2.0*sh - sr0
	sz0 = sv

	return sr0, s00, sz0

def yeild_DP(p, q, a, b):
	return q - b*p - a

def dFdpq_DP(p, q, param_b):
	dFdp = -param_b
	dFdq = 1.0
	dFda = -1.0

	return dFdp, dFdq, dFda

def dGdpq_DP(p, q, param_b_potential):
	dGdp = -param_b_potential
	dGdq = 1.0

	return dGdp, dGdq

def compute_q_ep_boudary_DP(sh,sv,param_a,param_b):
	p0,q0 = epwAnal.invariants(sh,sh,sv)
	return param_b*p0 + param_a

def compute_param_b(frictionAngle):
	sin_frictionAngle = np.sin(frictionAngle)
	return 6.0*sin_frictionAngle/(3.0-sin_frictionAngle)

def DruckerPragerModel(a0_a_ratio, sh, sv, nu, a0, G, defaultCohesion, defaultFrictionAngle, defaultDilationAngle, defaultHardeningRate):
	a = a0/a0_a_ratio
	xi_well = 1.0 - a0/a # the auxiliary variable xi at the wellbore, xi = (a-a0)/a = 1-1/(a/a0)
	nPoints = 1000

	defaultFrictionAngle *= np.pi/180.0 # converted to rad
	defaultDilationAngle *= np.pi/180.0 # converted to rad
	
	param_b = compute_param_b(defaultFrictionAngle)
	param_b_potential = compute_param_b(defaultDilationAngle)
	param_a = defaultCohesion * param_b / np.tan(defaultFrictionAngle)

	# ELastic trial
	pw = 2.0*G*xi_well + sh
	r_e,sr_e,s0_e,sz_e = epwAnal.solution_elastic(sh,sv,1.0,pw)
	p_e, q_e = epwAnal.invariants(sr_e[0],s0_e[0],sz_e[0])

	if (yeild_DP(p_e, q_e, param_a, param_b)<=0): # Pure elastic wellbore
		return [0],[0],[0],[0],r_e,sr_e,s0_e,sz_e,pw,p_e,q_e

	else: # Elastic-plastic wellbore
		q_ep_boundary = compute_q_ep_boudary_DP(sh,sv,param_a,param_b)
		# Elastic-Plastic boundary
		sr_ep_boundary, s0_ep_boundary, sz_ep_boundary = elastic_plastic_boundary_DP(sh,sv,q_ep_boundary)
		sr_p = [sr_ep_boundary]
		s0_p = [s0_ep_boundary]
		sz_p = [sz_ep_boundary]
		epsV_p = [0] # strain calculated starting from the reference state that is the ep boundary

		# Eq. 36 in Chen and Abousleiman 2017
		xi_ep_boundary = (sr_ep_boundary - sh) / 2.0 / G
	
		# Mesh from elastic-plastic boundary to wellbore
		dxi = (xi_well - xi_ep_boundary) / (nPoints - 1)
		xi = np.linspace(xi_ep_boundary, xi_well, nPoints)
	
		for i in range(1, nPoints):
			
			p,q = epwAnal.invariants(sr_p[i-1],s0_p[i-1],sz_p[i-1])
			
			exp_epsilon_V = np.exp(epsV_p[i-1])# epsV is the volumetric strain from the elastic-plastic boundary
			
			E = 2.0 * G * (1.0 + nu)

			param_a = q - param_b*p # update the hardening parameter with plastic yield condition F=0
			dFdp, dFdq, dFda = dFdpq_DP(p,q,param_b)
			dGdp, dGdq = dGdpq_DP(p,q,param_b_potential)
			dadLambda = defaultHardeningRate
			
			dFdHardningParam,dHardningParamdPlasticVar = dFda,dadLambda
			
			sr_i, s0_i, sz_i, epsV_i = epwAnal.solution_plastic(sr_p[i - 1], s0_p[i - 1], sz_p[i - 1], epsV_p[i-1], xi[i - 1], dxi, dGdp, dGdq, E, nu, dFdHardningParam,dHardningParamdPlasticVar,1.0)
			
			sr_p.append(sr_i)
			s0_p.append(s0_i)
			sz_p.append(sz_i)        
			epsV_p.append(epsV_i)

		# Wellbore surface stress
		pw = sr_i
		p = (sr_i + s0_i + sz_i)/3.0
		q = np.sqrt(0.5) * np.sqrt( (sr_i-s0_i)**2.0 + (sr_i-sz_i)**2.0 + (s0_i-sz_i)**2.0 )

		# Compute the normalized radial coordinate r/a
		r_p = epwAnal.compute_radialCoordinate(xi, epsV_p)
	
		# Elastic zone
		r_ep_boundary = r_p[0]
		r_e,sr_e,s0_e,sz_e = epwAnal.solution_elastic(sh,sv,r_ep_boundary,sr_ep_boundary)
			
		return r_p,sr_p,s0_p,sz_p,r_e,sr_e,s0_e,sz_e,pw,p,q


