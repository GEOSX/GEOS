import sys
import matplotlib
import numpy as np
import scipy.special

import math
import matplotlib.pyplot as plt


class HydrofractureSolutions:

    def epsilon_m(self, mup, Ep, t):
        raise RuntimeError('epsilon_m not defined ')
        return 0

    def Lm(self, Ep, Q0, t, mup):
        raise RuntimeError('Lm not defined ')
        return 0

    def epsilon_k(self, KIp, Ep, Q0, t):
        raise RuntimeError('epsilon_k not defined ')
        return 0

    def Lk(self, Q0, Ep, t, KIp):
        raise RuntimeError('Lk not defined ')
        return 0

    def gamma_m0(self):
        raise RuntimeError('gamma_m0 not defined ')
        return 0

    def gamma_k0(self):
        raise RuntimeError('gamma_k0 not defined ')
        return 0

    def Omega_m0_bar(self, rho):
        raise RuntimeError('Omega_m0_bar not defined ')
        return 0

    def PI_m0(self, rho):
        raise RuntimeError('PI_m0 not defined ')
        return 0

    def Omega_k0(self, rho):
        raise RuntimeError('Omega_k0 not defined ')
        return 0

    def PI_k0(self):
        raise RuntimeError('PI_k0 not defined ')
        return 0

    ## calcualte the assymtotic solutions as given by
    # Savitski and Detournay, International Journal of Solids and Structures (2002)
    def Solutions(self, mu, E, q, KI, time, time_radplots, distanceFromInjection):
        mup = 12 * mu

        PI = math.pi
        num_radplots = time_radplots.size

        # rho is the dimensionless radius with a range (0,1)
        numRho = 999
        rho = np.arange(numRho)
        rho = (rho + 1.0) / (numRho + 1.0)

        # Omega_m0 is the dimensionless crack opening for the zero-toughness solution
        # Omega_k0 is the dimensionless crack opening for the zero-viscosity solution
        # PI_m0 is the dimensinoless net pressure  for

        Omega_k0 = np.zeros(shape=[numRho])
        Omega_m0_bar = np.zeros(shape=[numRho])
        Omega_m0 = np.zeros(shape=[numRho])
        PI_m0 = np.zeros(shape=[numRho])
        PI_k0 = np.zeros(shape=[numRho])

        w_tough = np.zeros(shape=[num_radplots, numRho])
        p_tough = np.zeros(shape=[num_radplots, numRho])
        w_viscous = np.zeros(shape=[num_radplots, numRho])
        p_viscous = np.zeros(shape=[num_radplots, numRho])

        for j in range(0, len(rho)):
            Omega_m0_bar[j] = self.Omega_m0_bar(rho[j])

            # Equation (69)
            PI_m0[j] = self.PI_m0(rho[j])

            # Equation (84)
            PI_k0[j] = self.PI_k0()

            #Equation (85)
            Omega_k0[j] = self.Omega_k0(rho[j])

        # Equation 39
        Omega_m0 = Omega_m0_bar * self.gamma_m0()

        if 0:
            fig = plt.figure()
            plt.plot(rho, Omega_m0_bar)
            plt.plot(rho, Omega_k0)
            fig = plt.figure()
            plt.plot(rho, PI_m0)
            plt.plot(rho, PI_k0)
            plt.show()

        # plots for quantities vs radius
        for i in range(0, num_radplots):

            t = time_radplots[i]

            # Equation 15
            KIp = 4.0 * math.sqrt(2.0 / PI) * KI

            w_viscous[i] = self.epsilon_m(mup, E, t) * self.Lm(E, q, t, mup) * Omega_m0
            p_viscous[i] = self.epsilon_m(mup, E, t) * E * PI_m0

            w_tough[i] = self.epsilon_k(KIp, E, q, t) * self.Lk(q, E, t, KIp) * Omega_k0
            p_tough[i] = self.epsilon_k(KIp, E, q, t) * E * PI_k0

        # plots for opening values vs time
        PvsTime_Tdom = np.empty_like(time)
        RvsTime_Tdom = np.empty_like(time)
        w0vsTime_Tdom = np.empty_like(time)

        P_Vdom = np.empty_like(time)
        R_Vdom = np.empty_like(time)
        w0_Vdom = np.empty_like(time)

        for j in range(0, len(time)):
            t = time[j]
            if t == 0.0:
                t = 1.0e-10

            # Equation 27
            em = self.epsilon_m(mup, E, t)
            Lm = self.Lm(E, q, t, mup)

            # Equation 15
            KIp = 4.0 * math.sqrt(2.0 / PI) * KI

            # Equation 32
            ek = self.epsilon_k(KIp, E, q, t)
            Lk = self.Lk(q, E, t, KIp)

            # equation 14
            R_Vdom[j] = Lm * self.gamma_m0()
            RvsTime_Tdom[j] = Lk * self.gamma_k0()

            RrhoV = rho * R_Vdom[j]
            RrhoT = rho * RvsTime_Tdom[j]

            PI_k0_wb = np.interp(distanceFromInjection, RrhoT, PI_k0)
            PI_m0_wb = np.interp(distanceFromInjection, RrhoV, PI_m0)

            Omega_k0_wb = np.interp(distanceFromInjection, RrhoT, Omega_k0)
            Omega_m0_wb = np.interp(distanceFromInjection, RrhoV, Omega_m0)

            # Equation 12
            w0vsTime_Tdom[j] = ek * Lk * Omega_k0_wb
            w0_Vdom[j] = em * Lm * Omega_m0_wb

            #Equation 13
            P_Vdom[j] = em * E * PI_m0_wb
            PvsTime_Tdom[j] = ek * E * PI_k0_wb

        return [
            rho, w_tough, p_tough, w_viscous, p_viscous, PvsTime_Tdom, RvsTime_Tdom, w0vsTime_Tdom, P_Vdom, R_Vdom,
            w0_Vdom
        ]


class KGDSolutions(HydrofractureSolutions):

    def epsilon_m(self, mup, Ep, t):
        return pow(mup / (Ep * t), 1.0 / 3.0)

    def Lm(self, Ep, Q0, t, mup):
        return pow(Ep * pow(Q0, 3) * pow(t, 4) / mup, 1.0 / 6.0)

    def epsilon_k(self, KIp, Ep, Q0, t):
        return pow(pow(KIp / Ep, 4.0) / (Q0 * t), 1.0 / 3.0)

    def Lk(self, Q0, Ep, t, KIp):
        return pow(Ep * Q0 * t / KIp, 2.0 / 3.0)

    def gamma_m0(self):
        return 0.616

    def gamma_k0(self):
        return 2.0 / pow(math.pi, 2.0 / 3.0)

    def Omega_m0_bar(self, rho):
        A0 = math.sqrt(3.0)
        A1 = -0.156
        B1 = 0.0663
        return A0 * pow( 1.0 - rho * rho , 2.0/3.0 ) \
             + A1 * pow( 1.0 - rho * rho , 5.0/3.0 ) \
             + B1 * ( 4.0 * math.sqrt(1-rho*rho) + 2.0*rho*rho*math.log( (1.0-math.sqrt(1-rho*rho))/(1.0+math.sqrt(1-rho*rho)) ) )

    def PI_m0(self, rho):
        A1 = 1.61750
        A2 = 0.39650
        B = 0.06858
        PIStar1 = 1.0 / (3.0 * math.pi) * scipy.special.beta(0.5, 2.0 / 3.0) * scipy.special.hyp2f1(
            -1.0 / 6.0, 1, 0.5, rho * rho)
        PIStar2 = 3.0 / (20.0 * math.pi) * scipy.special.beta(
            -1.5, 8.0 / 3.0) * (2.0 / 3.0 * rho * rho * scipy.special.hyp2f1(-1.0 / 6.0, 2, 1.5, rho * rho) -
                                0.5 * scipy.special.hyp2f1(-7.0 / 6.0, 1, 0.5, rho * rho))

        return A1 * PIStar1 + A2 * PIStar2 + B * (2.0 - math.pi * abs(rho))


#        A0 =  math.sqrt(3.0)
#        A1 = -0.156
#        B1 =  0.0663
#        B = 2.58711
#        return  1.0 / ( 3.0 * math.pi ) * B * A0 * scipy.special.hyp2f1(-1.0/6.0,1,0.5,rho*rho) \
#              + 10.0/7.0 * A1 * scipy.special.hyp2f1(-7.0/6.0,1,0.5,rho*rho) \
#              + B1*(2-math.pi*abs(rho))

    def Omega_k0(self, rho):
        return self.gamma_k0() * pow(math.pi, 1.0 / 3.0) / 2.0 * math.sqrt(1 - pow(rho, 2))

    def PI_k0(self):
        return (pow(math.pi, 1.0 / 3.0) / 8.0)


class PennySolutions(HydrofractureSolutions):

    def epsilon_m(self, mup, Ep, t):
        return pow(mup / (Ep * t), 1.0 / 3.0)

    def Lm(self, Ep, Q0, t, mup):
        return pow(Ep * pow(Q0, 3) * pow(t, 4) / mup, 1.0 / 9.0)

    def epsilon_k(self, KIp, Ep, Q0, t):
        return pow(pow(KIp / Ep, 6.0) / (Q0 * t), 1.0 / 5.0)

    def Lk(self, Q0, Ep, t, KIp):
        return pow(Ep * Q0 * t / KIp, 2.0 / 5.0)

    def gamma_m0(self):
        return 0.6955

    def gamma_k0(self):
        return pow(3 / (math.sqrt(2.0) * math.pi), 2.0 / 5.0)

    def Omega_m0_bar(self, rho):
        B1 = 9.269e-2
        A1 = 3.581e-1
        C1 = 6.846e-1
        C2 = 7.098e-2
        ot1 = (math.sqrt(70.0) / 3.0 * C1 + 4.0 * math.sqrt(5.0) / 9.0 * C2 *
               (13 * rho - 6)) * pow(1.0 - rho, 2.0 / 3.0)
        ot2 = B1 * 8.0 / math.pi * (math.sqrt(1.0 - rho) - rho * math.acos(rho))
        return (ot1 + ot2)

    def PI_m0(self, rho):
        B1 = 9.269e-2
        A1 = 3.581e-1
        C1 = 6.846e-1
        C2 = 7.098e-2
        return A1 * (2.479 - 2.0 / 3.0 / pow(1.0 - rho, 1.0 / 3.0)) - B1 * (math.log(rho / 2.0) + 1)

    def Omega_k0(self, rho):
        return pow(3.0 / 8.0 / math.pi, 0.2) * math.sqrt(1 - pow(rho, 2))

    def PI_k0(self):
        return (math.pi / 8.0 * pow(math.pi / 12, 1.0 / 5.0))


class PKN_viscosityStorageDominated:

    def __init__(self, E, nu, KIC, mu, Q0, t, h):
        Ep = E / (1.0 - nu**2.0)
        self.t = t
        self.Q0 = Q0
        self.mu = mu
        self.Ep = Ep
        self.h = h

    def analyticalSolution(self):
        t = self.t
        Q0 = self.Q0
        mu = self.mu
        Ep = self.Ep
        h = self.h

        halfLength = 0.3817 * ((Ep * Q0**3.0 * t**4.0) / (mu * h**4.0))**(1.0 / 5.0)

        inletAperture = 3.0 * ((mu * Q0 * halfLength) / (Ep))**(1.0 / 4.0)

        inletPressure = ((16.0 * mu * Q0 * Ep**3.0 * halfLength) / (np.pi * h**4.0))**(1.0 / 4.0)

        return [halfLength, inletAperture, inletPressure]



if __name__ == "__main__":
    mu = 0.0005
    E = 3.2e10
    q = 0.02
    KI = 5.47e5

    queriesDone = 0

    time = np.arange(401.0)

    #hfsolns = PennySolutions()
    hfsolns = KGDSolutions()

    radTimes = np.array([100.0, time[-1]])
    print("radTimes = ", radTimes)

    solns = hfsolns.Solutions(mu, E, q, KI, time, radTimes)

    rho = solns[0]
    w_tough = solns[1]
    p_tough = solns[2]
    w_viscous = solns[3]
    p_viscous = solns[4]
    P_Tdom = solns[5]
    R_Tdom = solns[6]
    w0_Tdom = solns[7]
    P_Vdom = solns[8]
    R_Vdom = solns[9]
    w0_Vdom = solns[10]

    N1 = 4
    font = {'size': 10}

    matplotlib.rc('font', **font)

    fig1 = plt.figure(figsize=(6, 4))    # a new figure window
    fig1.subplots_adjust(wspace=0.4, bottom=0.1, left=0.11, right=0.6, top=0.95)

    ax1 = fig1.add_subplot(3, 1, 1)    # specify (nrows, ncols, axnum)
    plt.plot(time, P_Tdom / 1.0e6, 'b-', label='$\mu$ = 0')
    plt.plot(time, P_Vdom / 1.0e6, 'r-', label='$K_{Ic}$ = 0')
    #plt.plot(time[2::N1],aper0[2::N1]*1000,'ks',fillstyle='none', markersize=3, label='GEOS-Cartesian Mesh' )
    #plt.plot(time2[::N1],aper02[::N1]*1000,'ko',fillstyle='none', markersize=3, label='GEOS-Radial Mesh' )
    #plt.xlabel('time (s)')
    plt.ylabel('Pressure (MPa)')
    plt.ylim([0.0, 0.25])
    #plt.legend(loc='lower right')
    plt.xticks(np.arange(min(time), max(time), 100.0))
    plt.xlim([min(time), max(time)])

    ax2 = fig1.add_subplot(3, 1, 2)    # specify (nrows, ncols, axnum)
    plt.plot(time, w0_Tdom * 1000, 'b-', label='$\mu$ = 0')
    plt.plot(time, w0_Vdom * 1000, 'r-', label='$K_{Ic}$ = 0')
    #plt.plot(time[2::N1],aper0[2::N1]*1000,'ks',fillstyle='none', markersize=3, label='GEOS-Cartesian Mesh' )
    #plt.plot(time2[::N1],aper02[::N1]*1000,'ko',fillstyle='none', markersize=3, label='GEOS-Radial Mesh' )
    #plt.xlabel('time (s)')
    plt.ylabel('Aperture0 (mm)')
    #plt.legend(loc='lower right')
    plt.xticks(np.arange(min(time), max(time), 100.0))
    plt.xlim([min(time), max(time)])

    ax3 = fig1.add_subplot(3, 1, 3)    # specify (nrows, ncols, axnum)
    plt.plot(time, R_Tdom, 'b-', label='$\mu$ => 0')
    plt.plot(time, R_Vdom, 'r-', label='$K_{Ic}$ => 0')
    #plt.plot(time[2::N1],Radius[2::N1],'ks',fillstyle='none', markersize=3, label='GEOS-Cartesian Mesh' )
    #plt.plot(time2[::N1],Radius2[::N1],'ko',fillstyle='none', markersize=3, label='GEOS-Radial Mesh' )
    plt.xlabel('time (s)')
    plt.ylabel('Radius (m)')
    plt.xticks(np.arange(min(time), max(time), 100.0))
    plt.xlim([min(time), max(time)])

    plt.legend(loc='center right', bbox_to_anchor=(1.75, 1.1), prop={'size': 10})

    #plt.savefig(prefix + '_vsTime.eps')
    plt.show()

#=======================================================================================================================
# N1 = 8
# fig2 = plt.figure(figsize=(6, 4))  # a new figure window
# fig2.subplots_adjust(wspace=0.4,bottom=0.1,left=0.11,right=0.6,top=0.95)
#
# ax2 = fig2.add_subplot(2, 1, 1)  # specify (nrows, ncols, axnum)
# plt.plot(solns[0],w_tough*1000,  'b-', label='$\mu$ => 0' )
# plt.plot(solns[0],w_viscous*1000,'r-', label='$K_{Ic}$ => 0' )
# plt.plot(radA[2::N1]/np.amax(Radius),Aperture[2::N1]*1000,'ks',fillstyle='none', markersize=3,label='GEOS-Cartesian Mesh' )
# plt.plot(radA2/np.amax(Radius2),Aperture2*1000,'ko',fillstyle='none', markersize=3,label='GEOS-Radial Mesh' )
# #plt.xlabel('Radius (m)')
# plt.ylabel('Aperture (mm)')
# #plt.legend(loc='upper right')
# plt.xlim([0, max(radA/np.amax(Radius)) ])
#
#
# ax2 = fig2.add_subplot(2, 1, 2)  # specify (nrows, ncols, axnum)
# plt.plot(solns[0],p_tough/1.0e6,  'b-', label='$\mu$ = 0' )
# plt.plot(solns[0],p_viscous/1.0e6,'r-', label='$K_{Ic}$ = 0' )
# plt.plot(radP[2::N1]/np.amax(Radius),Pressure[2::N1]/1.0e6,'ks',fillstyle='none', markersize=3, label='GEOS-Cartesian Mesh' )
# plt.plot(radP2[2::N1]/np.amax(Radius2),Pressure2[2::N1]/1.0e6,'ko',fillstyle='none', markersize=3, label='GEOS-Radial Mesh' )
# plt.xlabel('Radial Coordinate')
# plt.ylabel('Pressure (MPa)')
# plt.xlim([0, max(radP/np.amax(Radius)) ])
# plt.legend(loc='center right', bbox_to_anchor=(1.75, 1.1),prop={'size':10})
#
#
#
# #ylim = plt.ylim()
# #pmin = np.amin(Pressure/1.0e6)
# #pmax = np.amax(Pressure/1.0e6)
# #range = np.ptp(Pressure/1.0e6)
# #plt.ylim([ pmin - range/5  , pmax + range/5 ])
# #plt.ylim([ -0.1  , 0.4 ])
#
# #plt.savefig(prefix + '_vsRadius.eps',dpi=400)
#
#
# plt.show()
#=======================================================================================================================
