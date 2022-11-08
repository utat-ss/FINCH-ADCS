from acstoolbox.constants.time_constants import *
from acstoolbox.foundation import so3

from numpy import pi as nppi
from numpy import mod as mod

# TODO:
# 1. Create a Cie(self, JD_TT_from_J2000)


class Earth:
    def __init__(self, clock):
        """Earth object describing terrestrial properties

        1. Earth-Fixed to Earth-Inertial Rotation (IAU2000)
              Entry Point: Cie
              Source: Astronomical Almanac
        """
        self.clock_ = clock

    def Cie(self, gregorian_utc):
        # Sidereal rotation.
        dut1 = self.clock_.GetdUTfromGregorian(gregorian_utc)
        gregorian_ut1 = gregorian_utc
        gregorian_ut1[5] += dut1
        jd_ut1 = self.clock_GregorianToJulianDate(gregorian_ut1)
        theta_rad = self.EarthRotationAngle(jd_ut1)

        return so3.C(3, -theta_rad)

    def EarthRotationAngle(self, jd_ut1):
        # Evaluate the Earth Rotation Angle for
        # Unit Test (Ephemerides-7): test_earth_rotation_angle
        djd_ut1 = jd_ut1 - JD_J2000
        theta_deg = mod(280.46061837504 + 360.985612288808 * djd_ut1, 360.0)

        return theta_deg / 180.0 * nppi
