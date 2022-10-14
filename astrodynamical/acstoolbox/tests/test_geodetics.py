""" test ellipsoid.py"""

# ACSToolbox packages.
from acstoolbox.geodetics import ellipsoid

# ACSToolbox packages.
import numpy as np
import pytest as pytest


# 1. Evaluate the site position vector at Mt.Evans, Colorado ( Geodetic latitude = 39.586667 degree, longitude = -105.640 degree, height = 4347.667 M
def test_ellipsoid():
    phi_d = 39.586667
    lamda = -105.640
    height = 4.347667

    # Calculate the site vector
    site_vector = ellipsoid.SitePositionVector2D(phi_d, height)

    # Example site vector
    site_vector_example = np.array([4925.42980258, 4045.49374256])

    # Test site vector in 2D
    assert site_vector == pytest.approx(site_vector_example, 1e-2)


# 2. Evaluate the site position vector at Toronto, Ontario ( Geodetic latitude = 43.6532  N, longitude = -79.2832 W, height = 0.076 M
def test_ellipsoid3D():
    phi_d = 43.6532
    lamda = -79.2832
    height = 0.076
    # Calculate the site vector
    site_vector = ellipsoid.SitePositionVector3D(phi_d, lamda, height)
    # Example Site Vector
    site_vector_example = np.array([859.52248524, -4541.5945528, 4380.34474178])

    # Test site vector in 3D
    assert site_vector == pytest.approx(site_vector_example, 1e-2)
