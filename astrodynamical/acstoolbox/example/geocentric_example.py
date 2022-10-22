""" ellipsoid .py example """
# ACSToolbox packages .
from acstoolbox . geodetics import ellipsoid

# python library .
import numpy as np
import pytest as pytest

# Evaluate the site position vector at Toronto , Ontario ( Geodetic latitude = 43.6532 N, longitude = -79.2832 W, height = 0.076 M)
if __name__ == " __main__ ":
    phi_d = 43.6532
    lamda = -79.2832
    height = 0.076
    # Calculate the site vector
    site_vector = ellipsoid . SitePositionVector3D ( phi_d , lamda , height )
    print (" Site Vector " = site_vector )
