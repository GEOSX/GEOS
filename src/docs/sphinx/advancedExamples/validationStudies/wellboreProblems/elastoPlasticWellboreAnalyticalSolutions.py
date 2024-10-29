import numpy as np

def invariants(sr,s0,sz):
    p = (sr + s0 + sz) / 3.0
    q = np.sqrt(((sr - s0)**2 + (sr - sz)**2 + (s0 - sz)**2) / 2.0)
    return p,q

def hardeningRate(dFdHardningParam,dHardningParamdPlasticVar,dPotentialdStressVar):
    return -dFdHardningParam*dHardningParamdPlasticVar*dPotentialdStressVar
    
def compute_a(sr,s0,sz,sigma,dGdp, dGdq):
    p,q = invariants(sr,s0,sz)
    dpdsigma = 1.0/3.0
    dqdsigma = 3.0*(sigma-p)/2.0/q
    return dGdp*dpdsigma + dGdq*dqdsigma

def compute_bij(sr, s0, sz, dGdp, dGdq, E, nu, dFdpc,dpcdepsVp,dPotentialdStressVar):
    h = hardeningRate(dFdpc,dpcdepsVp,dPotentialdStressVar)
    
    ar = compute_a(sr,s0,sz,sr,dGdp,dGdq)
    a0 = compute_a(sr,s0,sz,s0,dGdp,dGdq)
    az = compute_a(sr,s0,sz,sz,dGdp,dGdq)
    
    y = 1/(h+1e-10) # 1e-10 to avoid division to zero
    
    tmp1 = 1.0 / E**2.0
    b11 = tmp1 * (1.0 - nu**2.0 + E * a0 * a0 * y + 2.0 * E * nu * a0 * az * y + E * az * az * y)
    b12 = tmp1 * (-E * ar * (a0 + nu * az) * y + nu * (1.0 + nu - E * a0 * az * y + E * az * az * y))
    b13 = tmp1 * (-E * ar * (nu * a0 + az) * y + nu * (1.0 + nu + E * a0 * a0 * y - E * a0 * az * y))
    b22 = tmp1 * (1.0 - nu**2.0 + E * ar * ar * y + 2.0 * E * nu * ar * az * y + E * az * az * y)
    b23 = tmp1 * (nu + nu**2.0 + E * nu * ar * ar * y - E * a0 * az * y - E * nu * ar * (a0 + az) * y)
    b33 = tmp1 * (1.0 - nu**2.0 + E * ar * ar * y + 2.0 * E * nu * ar * a0 * y + E * a0 * a0 * y)
    b21 = b12
    b31 = b13
    b32 = b23

    delta = -(1.0 + nu) / E**3.0 * (-1.0 + nu + 2.0 * nu**2.0 + E * (-1.0 + nu) * ar * ar * y + E *
                                               (-1.0 + nu) * a0 * a0 * y - 2.0 * E * nu * a0 * az * y -
                                               E * az * az * y + E * nu * az * az * y - 2.0 * E * nu * ar *
                                               (a0 + az) * y)

    return b11,b12,b13,b22,b23,b33,b21,b31,b32,delta
    

def solution_plastic(sr, s0, sz, epsV, xd, dx, dGdp, dGdq, E, nu, dFdpc,dpcdepsVp,dPotentialdStressVar):
    b11,b12,b13,b22,b23,b33,b21,b31,b32,delta = compute_bij(sr, s0, sz, dGdp, dGdq, E, nu, dFdpc,dpcdepsVp,dPotentialdStressVar)

    exp_epsilon_V = np.exp(epsV) # epsV is the volumetric strain from the elastic-plastic boundary
    tmp = 1.0 - 2.0*xd - exp_epsilon_V # Use natural logarithm form for large volume strain
    
    DrDx = -(sr - s0) / tmp
    sr += DrDx * dx

    D0Dx = -b21 / b11 * ((sr - s0) / tmp + (b11 - b12) / delta /
                         (1 - xd)) - (b22 - b21) / delta / (1 - xd)
    s0 += D0Dx * dx

    DzDx = -b31 / b11 * ((sr - s0) / tmp + (b11 - b12) / delta /
                         (1 - xd)) - (b32 - b31) / delta / (1 - xd)
    sz += DzDx * dx

    DEpsDx = -delta / b11 * ((sr - s0) / tmp + (b11 - b12) / delta / (1 - xd))
    epsV += DEpsDx*dx
    
    return sr, s0, sz, epsV

def compute_radialCoordinate(xi,epsV):
    # Note that xi and r start from elastic-plastic boundary to the wellbore
    integral = 0.0 #integral from xi_well to xi_ep_boundary
    
    nPoints = len(xi)
    r = np.zeros(nPoints)
    
    # Loop from wellbore to the elastic-plastic boundary
    r[nPoints-1] = 1.0 # this is the normalized radial coordinate r/a
    for i in range(nPoints-1, 0, -1):
        dxi = xi[i-1]-xi[i]
        exp_epsilon_V = np.exp(epsV[i-1])
    
        integral += dxi / (1.0 - 2.0*xi[i - 1] - exp_epsilon_V) # Use natural logarithm form for large volume strain
        r[i - 1] = np.exp(integral)
        
    return r

def solution_elastic(sh,sv,r_ep_boundary,sr_ep_boundary):
    r_e = np.linspace(r_ep_boundary, 100, 1000)
    sr_e = sh + (sr_ep_boundary - sh) * (r_ep_boundary / r_e)**2.0
    s0_e = sh - (sr_ep_boundary - sh) * (r_ep_boundary / r_e)**2.0
    sz_e = sv*np.ones(len(r_e))
    
    return r_e,sr_e,s0_e,sz_e
