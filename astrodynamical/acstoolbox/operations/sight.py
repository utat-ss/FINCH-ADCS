""" This module contains functions to solve whether or not the line of sight (LOS) exists between two position vectors for the SIGHT and LIGHT problem.
    SIGHT: whether or not the LOS exists between two satellite vectors.
    LIGHT: whether or not the LOS exsits between the sun vector and satellite vector"""

# Standard packages.
import numpy as np

# ACSToolbox packages.
from acstoolbox.ephemerides.sun import Sun
from acstoolbox.time.clock import Clock


def VectortoScalar(r):
    value = 0
    length = r.size
    for i in range(length):
        value = value + (r[i] ** 2)
    scalar = np.sqrt(value)
    return scalar


def min_parametric(r1, r2):
    # find the minimum parametric value.
    s_r1 = VectortoScalar(r1)
    s_r2 = VectortoScalar(r2)
    tau_min = ((s_r1**2) - np.dot(r1, r2)) / (
        (s_r1**2) + (s_r2**2) - 2 * np.dot(r1, r2)
    )
    return tau_min


def sight(r1, r2):

    s_r1 = VectortoScalar(r1)
    tau_min = min_parametric(r1, r2)

    # Check
    if (tau_min < 0) or (tau_min > 1):
        return True
    else:
        c_vec = (1 - tau_min) * (s_r1**2) + np.dot(r1, r2) * tau_min

        if c_vec >= 1.0:
            return True
        else:
            return False


def light(r1, c_utc):

    clock = Clock()
    sun = Sun(clock)

    # find sun vector for sight
    s_mod = sun.GetUnitMODPositionFromUTC(c_utc) * 149597870.691
    LOS = sight(s_mod, r1)
    return LOS
