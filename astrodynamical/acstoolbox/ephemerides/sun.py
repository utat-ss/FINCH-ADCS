import numpy as np
from acstoolbox.constants.celestial_constants import *
from acstoolbox.constants.time_constants import *
from acstoolbox.ephemerides import astrodynamics


def Wrap(deg):
    return np.mod(deg, 360.0)


def WrapToRad(deg):
    return Wrap(deg) / 180.0 * np.pi


def DegToRad(deg):
    return deg / 180.0 * np.pi


# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
class Sun:
    def __init__(self, time):
        """Low Precision Sun Vector
        Source: Astronomical Almanac 1992:C24
        Accuracy: 0.01 deg
        Inputs: 1. Clock object, UTC time
        Output: 1. Sun vector in MOD frame
        Comments: This is a low precision algorithm with an accuracy of 0.01 degrees
                  in the Mean of Date (MOD) frame.
        """

        self.time_ = time

    def SunMeanLongitude(self, c_tdb):
        lambda_mean_deg = Wrap(280.460 + 36000.771 * c_tdb)
        return lambda_mean_deg

    def SunMeanAnomaly(self, c_tdb):
        theta_mean_rad = WrapToRad(357.5277233 + 35999.05034 * c_tdb)
        return theta_mean_rad

    def SunEclipticLatitude(self, c_tdb):
        lambda_mean_deg = self.SunMeanLongitude(c_tdb)
        theta_mean_rad = self.SunMeanAnomaly(c_tdb)
        lambda_ecliptic_rad = WrapToRad(
            lambda_mean_deg
            + 1.914666471 * np.sin(theta_mean_rad)
            + 0.019994643 * np.sin(2 * theta_mean_rad)
        )
        return lambda_ecliptic_rad  # in radians

    def GetUnitMODPositionFromTUT1(self, c_tdb):
        # An appropriate approximation to c_tdb is t_ut1.
        # Unit Test: 1. test_sun_vector_from_utc (indirectly)
        S_MOD = self.GeTKmMODposittionfromUTC(c_tdb)
        return S_MOD / np.linalg.norm(S_MOD)

    def GetMODFromJSJ2000UTC(self, js_j2000_utc):
        jd_utc = self.time_.JSJ2000ToJD(js_j2000_utc)
        mjd_utc = self.time_.MJD(jd_utc)
        dut1_s = self.time_.GetEarthObservationParameter(mjd_utc, "dUT1_s")
        t_ut1 = self.time_.JDToT(jd_utc + dut1_s / DAY_IN_SECONDS)

        return self.GetUnitMODPositionFromTUT1(t_ut1)

    # -----------------------------------------------------------------------
    def GetUnitMODPositionFromUTC(self, gregorian_utc):
        """Example Call:
        clock = acstoolbox.time.clock.Clock()
        sun = Sun(clock)
        s_mod = sun.GetMODPositionFromUTC(2006, 4, 2, 0, 0, 0)
        """

        # Unit Test: 1. test_sun_vector_from_utc

        t_ut1 = self.time_.UTCGregotianToTUT1(gregorian_utc)

        return self.GetUnitMODPositionFromTUT1(t_ut1)

    def GeTKmMODposittionfromUTC(self, c_tdb):
        lambda_mean_deg = self.SunMeanLongitude(c_tdb) # clock

        lambda_mean_rad = DegToRad(lambda_mean_deg)
        theta_mean_rad = self.SunMeanAnomaly(c_tdb)
        lambda_ecliptic_rad = self.SunEclipticLatitude(c_tdb)
        obliquity = WrapToRad(23.439291 - 0.0130042 * c_tdb)
        s_MOD = np.array(
            [
                np.cos(lambda_ecliptic_rad),
                np.cos(obliquity) * np.sin(lambda_ecliptic_rad),
                np.sin(obliquity) * np.sin(lambda_ecliptic_rad),
            ]
        )

        return s_MOD
