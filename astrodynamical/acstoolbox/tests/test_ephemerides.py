import sys
sys.path.append("C:\\Users\\khang\\OneDrive\\Music\\Documents\\ADCS\\simulation")
# ACSToolbox packages.
from acstoolbox.time.clock import Clock
from acstoolbox.ephemerides.sun import Sun
from acstoolbox.ephemerides.moon import Moon
from acstoolbox.ephemerides.earth import Earth

# Standard packages.
import numpy as np
import pytest as pytest

# TODO
# 1. Evaluate unit sun vector from JS from J2000 (in TT, not UTC)
# 2. Updated test_sun_vector_from_utc using current Astronomical Almanac

# 1. Evaluate unit sun vector from UTC Gregorian date.
def test_sun_vector_from_utc():
    clock = Clock()
    sun = Sun(clock)

    # Evaluate the sun vector in the MOD frame on 2 April 2006 00hh00mm00ss.
    epoch_gregorian_utc = [2006, 4, 2, 0, 0, 0]
    s_mod = sun.GetUnitMODPositionFromUTC(epoch_gregorian_utc)
    print(s_mod)

    # Astronomical Almanac (ICRS).
    s_mod_almanac = np.array([0.9776782, 0.1911521, 0.0828717])

    # Test the accuracy of each unit vector component.
    assert s_mod == pytest.approx(s_mod_almanac, 1e-2)

    # Test the angular difference of the unit sun vector.
    us_mod_almanac = s_mod_almanac / np.linalg.norm(s_mod_almanac)
    dphi = np.arccos(np.dot(s_mod, us_mod_almanac)) / np.pi * 180.0
    assert dphi < 1e-1


# 2. Evaluate moon vector from Julian Second from the J2000 epoch in ER
def test_moon_vector_from_tt():

    clock = Clock()
    moon = Moon(clock)
    # Evaluate the moon vector from April 28, 1994, 0:00 UTC
    gregorian_utc = [1994, 4, 28, 0, 0, 0]
    jsj2000_tt = clock.GregorianToJSJ2000(gregorian_utc)
    m_vector = moon.VectorFromTT(jsj2000_tt)

    # Astronomical Almanac (ICRS).
    m_mod_almanac = np.array([-21.0469963, -48.84993696, -19.86376033])

    assert m_vector == pytest.approx(m_mod_almanac, 1e-2)


# 3. Evaluate moon vector from UTC in ER
def test_moon_vector_from_utc():

    clock = Clock()
    moon = Moon(clock)

    # Evaluate the moon vector from April 28, 1994, 0:00 UTC
    gregorian_utc = [1994, 4, 28, 0, 0, 0]

    m_vector = moon.VectorFromUTC(gregorian_utc)

    # Astronomical Almanac (ICRS).
    m_mod_almanac = np.array([-21.0469963, -48.84993696, -19.86376033])

    assert m_vector == pytest.approx(m_mod_almanac, 1e-2)


# 4. Evaluate moon vector from UTC in KM
def test_moon_vector_in_km():
    clock = Clock()
    moon = Moon(clock)

    # Evaluate the moon vector from April 28, 1994, 0:00 UTC
    gregorian_utc = [1994, 4, 28, 0, 0, 0]
    m_vector = moon.VectorFromUTCinKM(gregorian_utc)

    # Example Value:
    example_vec = [-134241.192, -311571.349, -126693.681]
    assert m_vector == pytest.approx(example_vec, 1e-2)


# 5. Evaluate normalize moon vector from UTC
def test_moon_vec_normal():
    clock = Clock()
    moon = Moon(clock)
    # Evaluate the moon vector from April 28, 1994, 0:00 UTC
    gregorian_utc = [1994, 4, 28, 0, 0, 0]
    m_vector = moon.VectorFromUTC(gregorian_utc)
    moon_vector = m_vector / np.linalg.norm(m_vector)
    # Check with the new fuction
    moonvector = moon.MoonVectorNormal(gregorian_utc)
    assert moon_vector == pytest.approx(moonvector, 1e-2)


# 6. Evaluate moon phase from April 3, 2020 at 0:00 UTC
def test_moon_phase():
    clock = Clock()
    moon = Moon(clock)
    # Evaluate the moon phase on April 3, 2020 at 0:00 UTC
    gregorian_utc = [2020, 4, 3, 0, 0, 0]
    phase = moon.PercentageMoonIlluminated(gregorian_utc)
    # Test moon phase
    phase_test = 66.7277
    assert phase == pytest.approx(phase_test, 1e-2)

# 7. Evaluate Earth rotation angle.
# Source: Astronomical Almanac 2018, pp. B23
def test_earth_rotation_angle():
    clock = Clock()
    earth = Earth(clock)

    # Evaluate earth rotation angle for September 1, 2018, 0:00 UTC
    gregorian_ut1 = [2018, 9, 1, 0, 0, 0]

    jd_ut1 = clock.GregorianToJulianDate(gregorian_ut1)
    tabulated_jd_ut1 = 2458362.5
    assert jd_ut1 == tabulated_jd_ut1

    # Earth rotation angle is accurate to 1e-5 rad (2 arcseconds)
    theta_rad = earth.EarthRotationAngle(jd_ut1)
    tabulated_theta_rad = np.deg2rad(339 + 52 / 60 + 20.6127 / 3600)
    assert theta_rad == pytest.approx(tabulated_theta_rad, 1e-5)


if __name__ == "__main__":
    test_sun_vector_from_utc()
    epoch_gregorian_utc = [1994, 4, 2, 0, 0, 0]
    clock = Clock()
    sun = Sun(clock)
    c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
    sun_vector_test = sun.GeTKmMODposittionfromUTC(c_tdb)
    sun_vector = sun_vector_test * 149597870.691
    print(sun_vector)
