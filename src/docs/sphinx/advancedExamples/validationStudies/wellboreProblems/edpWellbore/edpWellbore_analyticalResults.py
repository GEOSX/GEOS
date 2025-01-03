import sys
sys.path.append('../')
import numpy as np
import elastoPlasticWellboreAnalyticalSolutions as epwAnal

#Note: this solution is developed by Chen and Abousleiman (2017) for the case with zero cohesion

def elastic_plastic_boundary_EDP(sh,sv,q_ep_boundary):
	# See eq. 25 of Chen and Abousleiman 2013
	sr0 = sh - np.sqrt(sh**2 - (4.0 * sh**2 + sv**2 - 2 * sh * sv - q_ep_boundary**2) / 3.0)
	s00 = 2.0*sh - sr0
	sz0 = sv

	return sr0, s00, sz0

def dFdpq_EDP(p, q):
	#F = q-tan(beta)*p for zero cohesion
	#F=0 -> tan(beta) = q/p
	param_b = q/p
	cos2Beta = 1.0/(param_b**2.0 + 1.0)

	dFdp = -param_b
	dFdq = 1.0
	dFdbeta = -p/cos2Beta

	return dFdp, dFdq, dFdbeta

def hardening_EDP(param_b, param_m, param_b_f, param_b_i):
	# gammaP is the deviatoric plastic strain
	# param_b = gammaP/(c + gammaP) * (param_b_f - param_b_i) + param_b_i
	# that yields: gammaP/(c + gammaP) = (param_b-param_b_i)/(param_b_f - param_b_i)
	# or c/(c + gammaP) = (param_b_f-param_b)/(param_b_f - param_b_i)
	# The derivative of the hardening law gives: dparam_bdGammaP = (param_b_f-param_b_i)c/(c + gammaP)**2
	# or dparam_bdGammaP = (param_b_f-param_b)**2/(param_b_f - param_b_i)/c
	# Finally dBetadGammaP = cos2Beta*(param_b_f-param_b)**2/(param_b_f - param_b_i)/c

	# Validation against Chen and Abousleiman 2017
	# h = -dFdbeta*dBetadGammaP
	# h = p*(param_b_f-param_b)**2/(param_b_f - param_b_i)/c
	# y = 1/h = c(param_b_f - param_b_i)p/(p*param_b_f - q)**2.0 that is identical to Chen and Abousleiman 2017

	cos2Beta = 1.0/(param_b**2+1.0)
	return cos2Beta*(param_b_f-param_b)**2/(param_b_f - param_b_i)/param_m

def compute_q_ep_boudary_EDP(sh,sv,param_b_i):
	p0,q0 = epwAnal.invariants(sh,sh,sv)
	return param_b_i*p0

def compute_param_b(frictionAngle):
	sin_frictionAngle = np.sin(frictionAngle)
	return 6.0*sin_frictionAngle/(3.0-sin_frictionAngle)

def EDP(a0_a_ratio, sh, sv, nu, a0, G, initialFrictionAngle, finalFrictionAngle, param_m):
	a = a0/a0_a_ratio
	xi_well = 1.0 - a0/a # the auxiliary variable xi at the wellbore, xi = (a-a0)/a = 1-1/(a/a0)
	nPoints = 1000

	initialFrictionAngle *= np.pi/180.0 # converted to rad
	finalFrictionAngle *= np.pi/180.0 # converted to rad

	param_b_i = compute_param_b(initialFrictionAngle)
	param_b_f = compute_param_b(finalFrictionAngle)

	# ELastic trial
	pw = 2.0*G*xi_well + sh
	r_e,sr_e,s0_e,sz_e = epwAnal.solution_elastic(sh,sv,1.0,pw)
	p_e, q_e = epwAnal.invariants(sr_e[0],s0_e[0],sz_e[0])

	if ((q_e/p_e)<param_b_i): # Pure elastic wellbore
		return [0],[0],[0],[0],r_e,sr_e,s0_e,sz_e,pw,p_e,q_e

	else: # Elastic-plastic wellbore
		q_ep_boundary = compute_q_ep_boudary_EDP(sh,sv,param_b_i)
	
		# Elastic-Plastic boundary
		sr_ep_boundary, s0_ep_boundary, sz_ep_boundary = elastic_plastic_boundary_EDP(sh,sv,q_ep_boundary)
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

			param_b = q/p # yield condition F=0
			dFdp, dFdq, dFdbeta = dFdpq_EDP(p,q)
			dBetadGammaP = hardening_EDP(param_b,param_m,param_b_f,param_b_i)
			dGdq = dFdq # assicated model
			dFdHardningParam,dHardningParamdPlasticVar,dPotentialdStressVar = dFdbeta,dBetadGammaP,dGdq
			
			sr_i, s0_i, sz_i, epsV_i = epwAnal.solution_plastic(sr_p[i - 1], s0_p[i - 1], sz_p[i - 1], epsV_p[i-1], xi[i - 1], dxi, dFdp, dFdq, E, nu, dFdHardningParam,dHardningParamdPlasticVar,dPotentialdStressVar)
			
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


