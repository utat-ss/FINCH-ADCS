""" Moon Position Vector"""

# Standard packages.
import numpy as np

# ACSToolbox packages.
from acstoolbox.ephemerides import astrodynamics
from acstoolbox.constants.celestial_constants import *
from acstoolbox.constants.time_constants import *
from acstoolbox.time.clock import Clock
from acstoolbox.ephemerides.sun import Sun


class Moon:
    def __init__(self, clock):
        """Low Precision Moon Vector
        Source:
        Accuracy:

        Initialization: Clock Object

        Input: UTC time, Terrestrial time (JS from J2000)
        Output: Moon Vector in [TBD] Frame
        """
        self.clock_ = clock

    # Unit Test (test_ephemerides.py): 2. test_moon_vector_from_tt

    def JulianSecondsfromJ2000toT_TBD(self, jsj2000_tt):
        # Initialize Clock.
        JD_J2000 = self.clock_.JSJ2000ToJD(jsj2000_tt)
        c_tdb = self.clock_.JDToT(JD_J2000)
        return c_tdb

    def Gregoriantoc_tdb(self, GregorianUT1):
        c_tdb = self.clock_.UTCGregotianToTUT1(GregorianUT1)
        return c_tdb

    def HorizontalParallax(self, c_tdb):
        # Calculating moon Horizontal Parallax
        mean_anomalies = astrodynamics.MoonMeanAnomaliesRad(c_tdb)
        mean_argument_Latitude = astrodynamics.MoonMeanArgumentofLatitudeRad(c_tdb)
        mean_elongation_from_sun = astrodynamics.MeanElongationfromSun(c_tdb)
        h_parallax = (
            0.9508
            + 0.0518 * np.cos(mean_anomalies)
            + 0.0095 * np.cos(mean_anomalies - 2 * mean_elongation_from_sun)
            + 0.0078 * np.cos(2 * mean_elongation_from_sun)
            + 0.0028 * np.cos(2 * mean_anomalies)
        )
        return h_parallax

    def EclipticLongitude(self, c_tdb):
        m_longitude = astrodynamics.MoonLongitude(c_tdb)  # Find Moon Longitude
        # Moon Variable
        mean_anomalies = astrodynamics.MoonMeanAnomaliesRad(c_tdb)
        mean_argument_Latitude = astrodynamics.MoonMeanArgumentofLatitudeRad(c_tdb)
        mean_elongation_from_sun = astrodynamics.MeanElongationfromSun(c_tdb)
        mean_anomaly_sun = astrodynamics.SunMeanAnomalies(c_tdb)
        # Calculate Ecliptic Longitude
        e_longitude = (
            m_longitude
            + 6.29 * np.sin(mean_anomalies)
            - 1.27 * np.sin(mean_anomalies - 2 * mean_elongation_from_sun)
            + 0.66 * np.sin(2 * mean_elongation_from_sun)
            + 0.21 * np.sin(2 * mean_anomalies)
            - 0.19 * np.sin(mean_anomaly_sun)
            - 0.11 * np.sin(2 * mean_argument_Latitude)
        )

        # If we get a negative Ecliptic Longitude, find the same location with positive angle.
        if e_longitude < 0:
            e_longitude = 360 - e_longitude
        return e_longitude

    def MagnitudeofPositionVector(self, c_tdb):
        # We can find the Magnitude of the Moon Position Vector using the Horizontal Parallex
        h_parallax = self.HorizontalParallax(c_tdb)
        mag_pos_vec = 1 / (
            np.sin(np.radians(h_parallax))
        )  # Magnitude of Moon Position Vector = 1 / sin (Horizontal Parallex)
        return mag_pos_vec

    def ObliquityofEcliptic(self, c_tdb):
        epsilon = 23.439291 - 0.0130042 * c_tdb
        return epsilon

    def EclipticLatitude(self, c_tdb):
        mean_anomalies = astrodynamics.MoonMeanAnomaliesRad(c_tdb)
        mean_argument_Latitude = astrodynamics.MoonMeanArgumentofLatitudeRad(c_tdb)
        mean_elongation_from_sun = astrodynamics.MeanElongationfromSun(c_tdb)

        # Calculating Ecliptic Latitude.
        e_latitude = (
            5.13 * np.sin(mean_argument_Latitude)
            + 0.28 * np.sin(mean_anomalies + mean_argument_Latitude)
            - 0.28 * np.sin(mean_argument_Latitude - mean_anomalies)
            - 0.17 * np.sin(mean_argument_Latitude - 2 * mean_elongation_from_sun)
        )
        return e_latitude

    def MoonVector(self, c_tdb):

        # Calculating ecliptic latitude, longtitude,parallax, obliquity and magnitude
        lamda = np.radians(self.EclipticLongitude(c_tdb))
        phi = np.radians(self.EclipticLatitude(c_tdb))
        mag = self.MagnitudeofPositionVector(c_tdb)
        epsilon = np.radians(self.ObliquityofEcliptic(c_tdb))

        # Moon Vector
        r_x = mag * (np.cos(phi) * np.cos(lamda))
        r_y = mag * (
            np.cos(epsilon) * np.cos(phi) * np.sin(lamda)
            - np.sin(epsilon) * np.sin(phi)
        )
        r_z = mag * (
            np.sin(epsilon) * np.cos(phi) * np.sin(lamda)
            + np.cos(epsilon) * np.sin(phi)
        )

        # Combine everything
        Moon_Vector = np.array([r_x, r_y, r_z])

        return Moon_Vector

    # Unit Test (test_ephemerides.py): 3. test_moon_vector_from_utc
    def VectorFromUTC(self, gregorian_utc):
        # Moon vector in ER
        c_tdb = self.Gregoriantoc_tdb(gregorian_utc)
        return self.MoonVector(c_tdb)

    def VectorFromTT(self, jsj2000_tt):
        # Moon vector in ER
        c_tdb = self.JulianSecondsfromJ2000toT_TBD(jsj2000_tt)
        return self.MoonVector(c_tdb)

    def VectorFromUTCinKM(self, gregorian_utc):
        # Moon vector in KM
        MoonVecinKm = self.VectorFromUTC(gregorian_utc) * 6378.14
        return MoonVecinKm

    def MoonVectorNormal(self, gregorian_utc):
        moonvector = self.VectorFromUTC(gregorian_utc)
        return moonvector / np.linalg.norm(moonvector)

    def MoonPhase(self, gregorian_utc):
        c_tdb = self.Gregoriantoc_tdb(gregorian_utc)
        sun = Sun(self.clock_)
        Lamda_sun = sun.SunEclipticLatitude(c_tdb)
        Lamda_moon = np.radians(self.EclipticLongitude(c_tdb))
        phase = Lamda_sun - Lamda_moon
        return phase

    def PercentageMoonIlluminated(self, gregorian_utc):
        MoonPhase = self.MoonPhase(gregorian_utc)
        PercentDisk = (100 / 2) * (1 - np.cos(MoonPhase))
        return PercentDisk
