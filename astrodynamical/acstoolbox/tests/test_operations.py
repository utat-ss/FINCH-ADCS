# Standard packages.
import numpy as np
import pytest

# ACSToolbox packages.
from acstoolbox.operations import sight

#TODO
# 1. Evaluate the existence of the LOS between two satellite vectors.
# 2. Evaluate the existence of the LOS between the sun vector and the satellite vector.


# 1. Evaluate there is no LOS between these two satellite position vectors.
def test_sight():
    r1 = np.array([0, -4464.696, -5102.509])
    r2 = np.array([0, 5740.323, 3189.068])
    # Conver to canonical units
    r1 = r1 / np.linalg.norm(r1)
    r2 = r2 / np.linalg.norm(r2)
    #test the sight algorithm
    LOS = sight.sight(r1, r2)
    assert LOS is False

# 2. Evaluate there is LOS between the Sun vector and this position vector.
def test_light():
    # Evaluate the sun vector in the MOD frame on 15 Febuary 1995 12hh00mm00ss.
    c_utc = [1995, 2, 15, 12, 0, 0]
    r1 = np.array([0, -4464.696, -5102.509])
    #test the light algorithm
    LOS = sight.light(r1, c_utc)
    assert LOS is True
    
