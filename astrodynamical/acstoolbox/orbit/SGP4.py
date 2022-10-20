import math
import numpy as np
from acstoolbox.constants.celestial_constants import *
from acstoolbox.constants.time_constants import *
from acstoolbox.time import time
from acstoolbox.constants.time_constants import *

def strsign2float(str):
    if str == '-':
        return -1.0

    # '+' or ' '
    return +1.0

def TLEParametersFromFile(tle_filepath, century):
    tle_param = {}

    with open(tle_filepath) as opened_file:
        tle_contents = opened_file.readlines()

        tle_param['bstar'] = strsign2float(tle_contents[0][53])*float("0."+tle_contents[0][54:59]) * 10**(strsign2float(tle_contents[0][59])*float(tle_contents[0][60]))
        tle_param['inclination'] = float(tle_contents[1][8:15])/180.0*np.pi
        tle_param['argument_perigee'] = float(tle_contents[1][34:42])/180.0*np.pi
        tle_param['eccentricity'] = float("0."+tle_contents[1][26:34])
        tle_param['right_ascension'] = float(tle_contents[1][17:25])/180.0*np.pi
        tle_param['mean_anomaly'] = float(tle_contents[1][43:51])/180.0*np.pi
        tle_param['mean_motion'] = float(tle_contents[1][52:63])*2*np.pi/1440
        tle_param['year'] = float(str(int(century/100)) + tle_contents[0][18:20])
        tle_param['fractional_days'] = float(tle_contents[0][20:32])

        jd_utc = time.YearFractionalDaystoJDUTC(tle_param['year'], tle_param['fractional_days'])
        tle_param['epoch_jsj2000_utc'] = (jd_utc - kJDJ2000) * kDayInSeconds



    return tle_param


class SGP4():
    def __init__(self, tle_param):
        """ SGP4 Algorithm
        Accuracy: 10s of km
        Source: http://www.celestrak.com/NORAD/documentation/spacetrk.pdf
                http://www.celestrak.com/publications/aiaa/2006-6753/AIAA-2006-6753.pdf
        Input:  1. TLE parameters converted from file
                   a) bstar
                   b) inclination [rad]
                   c) argument of perigee [rad]
                   d) eccentricity
                   e) RAAN [rad]
                   f) mean anomaly [rad]
                   g) mean motion [rad/min]
                2. Time elapsed in [min] from TLE epoch
        Output: 1. Position [km] in TEME frame
                2. Velocity [km/s] in TEME frame
        Example:
        # 1 88888U 80 275.98708465 .00073094 13844-3 66816-4 0 8
        # 2 88888 72.8435 115.9689 0086731 52.6988 110.5714 16.05824518 105
        # At 0 min:
        # r     2328.96594238 -5995.21600342 1719.97894287
        # v     2.91110113 -0.98164053 -7.09049922
        # At 360 min:
        # r     2456.00610352 -6071.94232177 1222.95977784
        # v     2.67852119 -0.44705850 -7.22800565
        """

        # Iterative solver for Kepler's equation.
        self.dEw0 = 1
        self.dEw_max = 1e-5
        self.n_iter_max = 30

        # TLE parameters read from file.
        self.bstar = tle_param['bstar']
        self.i0_rad = tle_param['inclination']
        self.w_rad = tle_param['argument_perigee']
        self.e0 = tle_param['eccentricity']
        self.Omega_rad = tle_param['right_ascension']
        self.M0_rad = tle_param['mean_anomaly']
        self.n0_radpmin = tle_param['mean_motion']

    def GetOrbitState(self, dt_min):
        """Sample Call
            tle_param = {}
            tle_param['bstar'] = 0.66816*(10**(-4))
            tle_param['inclination'] = 72.8435/180.0*np.pi
            tle_param['argument_perigee'] = 52.6988/180.0*np.pi
            tle_param['eccentricity'] = 0.0086731
            tle_param['right_ascension'] = 115.9689/180.0*np.pi
            tle_param['mean_anomaly'] = 110.5714/180.0*np.pi
            tle_param['mean_motion'] = 16.05824518*2*np.pi/1440
            sgp4 = tle.SGP4(tle_param)
            sgp4.GetOrbitState(0.0)
            sgp4.GetOrbitState(360.0)
        """
        # Astrodynamic constants based on Earth model.
        ke = (3600.0 * kG_E / (kr_E**3))**0.5
        k2 = 5.41308*(10**(-4))
        A30 = -kJ3
        k4 = 0.62098875*(10**(-6))
        Theta = np.cos(self.i0_rad)

        # Astrodynamical constants assuming perigee above 156km.
        q0 = 1.0 + 120/kr_E
        s = 1 + 78/kr_E

        # Recover the original mean motion (ni) and semi-major axis (ai).
        a1 = (ke/self.n0_radpmin) ** (2.0/3.0)
        d1 = 1.5 * k2/(a1**2) * (3*(Theta**2) - 1)/((1 - self.e0**2)**1.5)
        a0 = a1 * (1 - (1.0/3.0)*d1 - (d1**2) - (134.0/81.0)*(d1**3))
        d0 = 1.5 * k2/(a0**2) * (3*(Theta**2) - 1)/((1 - self.e0**2)**1.5)
        ni = self.n0_radpmin / (1 + d0)
        ai = a0 / (1 - d0)

        # Astrodynamical constants for evaluating secular effects of drag and gravitation.
        Tsi = 1 / (ai - s)
        B0 = (1-(self.e0**2))**0.5
        eta = ai * self.e0 * Tsi

        # Exponential constants to simplify secular expressions.
        k22 = k2**2
        kc2 = ai**2
        kc3 = ai**4
        kc4 = (q0-s)**4
        kc5 = Tsi**4
        kc6 = Theta**2
        kc7 = Theta**4
        kc8 = eta**2
        kc9 = eta**3
        kc10 = 1-(eta**2)
        kc11 = B0**2
        kc12 = B0**4
        kc13 = B0**8
        kc14 = kc10**(-3.5)

        # Coefficients used to evaluate the secular effects of drag and gravitation.
        C2 = kc4 * kc5 * ni * kc14 * (ai*(1+1.5*kc8 + 4*self.e0*eta
                                      + self.e0*(eta**3))
                                      + 1.5*(k2*Tsi/kc10)*(-0.5 + 1.5*kc6)*(8.0+24.0*kc8+3*(eta**4)))
        C1 = self.bstar * C2
        C3 = kc4*(Tsi**5)*A30*ni*kAE*np.sin(self.i0_rad)/(k2*self.e0)
        C4 = 2.0*ni*kc4*kc5*ai*kc11*kc14*((2*eta*(1+self.e0*eta)+0.5*self.e0+0.5*kc9)
                                           - 2.0*k2*Tsi/(ai*kc10)*(3.0*(1.0-3.0*kc6)*(1.0 + 1.5*kc8-2.0*self.e0*eta
                                           - 0.5*self.e0*kc9)+0.75*(1-kc6)*(2*kc8-self.e0*eta-self.e0*kc9)*np.cos(2*self.w_rad)))
        C5 = 2.0*kc4*kc5*ai*kc11*kc14 * \
            (1.0+(11.0/4.0)*eta*(eta+self.e0)+self.e0*kc9)
        D2 = 4.0*ai*Tsi*(C1**2)
        D3 = (4.0/3.0)*ai*(Tsi**2)*(17.0*ai + s)*(C1**3)
        D4 = (2.0/3.0)*ai*(Tsi**3)*(221.0*ai + 31.0*s)*(C1**4)

        # Secular effects of atmospheric drag and gravitation.
        MDF = self.M0_rad+(1.0+3.0*k2*(-1.0+3.0*kc6)/(2.0*kc2*(B0**3))+3.0*k22*(13.0-78.0*kc6
                           + 137.0*kc7)/(16.0*kc3*(B0**7)))*ni*dt_min
        wDF = self.w_rad+(-3.0*k2*(1.0-5.0*kc6)/(2.0*kc2*kc12)+3.0*k22*(7.0-114.0*kc6
                          + 395.0*kc7)/(16.0*kc3*kc13)+5.0*k4*(3.0-36.0*kc6+49.0*kc7)/(4.0*kc3*kc13))*ni*dt_min
        OmegaDF = self.Omega_rad+(-3.0*k2*Theta/(kc2*kc12)+3.0*k22*(4.0*Theta-19.0*(Theta**3))/(2.0*kc3*kc13)
                                  + 5.0*k4*Theta*(3.0-7.0*kc6)/(2.0*kc3*kc13))*ni*dt_min
        dw = self.bstar*C3*np.cos(self.w_rad)*dt_min
        dM = -(2.0/3.0)*kc4*self.bstar*kc5*(kAE/(self.e0*eta))*(((1+eta*np.cos(MDF))**3)-((1+eta*np.cos(self.M0_rad))**3))
        Mp = MDF + dw + dM
        w = wDF - dw - dM
        Omega = OmegaDF - (10.5)*(ni*k2*Theta/(kc2*kc11))*C1*(dt_min**2)
        e = self.e0 - self.bstar*C4*dt_min - self.bstar * \
            C5*(np.sin(Mp)-np.sin(self.M0_rad))
        a = ai*((1-C1*dt_min-D2*(dt_min**2)-D3*(dt_min**3)-D4*(dt_min**4))**2)
        IL = Mp + w + Omega + ni*(1.5*C1*(dt_min**2)+(D2+2*(C1**2))*(dt_min**3)
                                  + 0.25*(3.0*D3+12.0*C1*D2+10.0*(C1**3))*(dt_min**4)
                                  + 0.2*(3.0*D4+12.0*C1*D3+6.0*(D2**2) + 30.0*(C1**2)*D2+15.0*(C1**4))*(dt_min**5))
        B = (1-(e**2))**0.5
        n = ke / (a**1.5)

        # Long-period periodic terms.
        axN = e*np.cos(w)
        IL_L = (A30*np.sin(self.i0_rad)/(8.0*k2*a*(B**2)))*axN*((3.0+5.0*Theta)/(1.0+Theta))
        ayNL = A30*np.sin(self.i0_rad) / (4.0*k2*a*(B**2))
        IL_T = IL + IL_L
        ayN = e*np.sin(w) + ayNL

        # Set up and solve Kepler's equation for Ew_k = E + w.
        U = IL_T - Omega
        Ew_k = U
        dEw = self.dEw0
        n_iter = 0
        while (dEw > self.dEw_max and n_iter < self.n_iter_max):
            Ew_kp1 = Ew_k + (U - ayN*np.cos(Ew_k) + axN*np.sin(Ew_k)- Ew_k)/(-ayN*np.sin(Ew_k) - axN*np.cos(Ew_k) + 1.0)
            dEw = abs(Ew_kp1 - Ew_k)
            Ew_k = Ew_kp1
            n_iter = n_iter + 1

        # Preliminary quantities for short-period periodics
        ecosE = axN*np.cos(Ew_k) + ayN*np.sin(Ew_k)
        esinE = axN*np.sin(Ew_k) - ayN*np.cos(Ew_k)
        eL = (((axN**2) + (ayN**2)))**0.5
        pL = a * (1.0-(eL**2))
        r = a * (1.0-ecosE)
        rdot = ke*((a**0.5)/r)*esinE
        rfdot = ke*(pL**0.5)/r
        cosu = a/r*(np.cos(Ew_k) - axN + ayN*esinE/(1.0+(1.0-(eL**2))**0.5))
        sinu = a/r*(np.sin(Ew_k) - ayN - axN*esinE/(1.0+(1.0-(eL**2))**0.5))
        u = math.atan2(sinu, cosu)
        dr = (k2/(2.0*pL))*(1-kc6)*np.cos(2*u)
        du = -(k2/(4.0*(pL**2)))*(7.0*kc6-1)*np.sin(2*u)
        dOmega = (3.0*k2*Theta/(2.0*(pL**2)))*np.sin(2*u)
        di = (3.0*k2*Theta/(2.0*(pL**2)))*np.sin(self.i0_rad)*np.cos(2*u)
        drdot = -(k2*n/pL)*(1.0-kc6)*np.sin(2*u)
        drfdot = (k2*n/pL)*((1.0-kc6)*np.cos(2*u)-1.5*(1.0-3.0*kc6))

        # Sum the short-period periodic terms to obtain osculating quantities.
        rk = r*(1.0-1.5*k2*((1.0-(eL**2))**0.5/((pL**2)))
                * (3.0*(Theta**2)-1.0))+dr
        uk = u + du
        Omegak = Omega + dOmega
        ik = self.i0_rad + di
        rkdot = rdot + drdot
        rfkdot = rfdot + drfdot

        # Unit orientation vectors.
        M = np.array([-np.sin(Omegak)*np.cos(ik), np.cos(Omegak)*np.cos(ik), np.sin(ik)])
        N = np.array([np.cos(Omegak), np.sin(Omegak), 0.0])
        U = M*np.sin(uk) + N * np.cos(uk)
        V = M*np.cos(uk) - N * np.sin(uk)

        # Calculate the position and velocity vectors in TEME.
        # Transform the velocity from [km/min] to [km/s].
        kunnormalized = kr_E/kAE
        r_teme = rk * U * kunnormalized
        v_teme = (rkdot*U + rfkdot*V) * kunnormalized * \
                 (kDayInMinutes/kDayInSeconds)

        return (r_teme, v_teme)