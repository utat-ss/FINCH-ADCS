"""ACS Toolbox: SO3 Module
This module contains utility functions for the SO(3) 3D Rotational Group.
These special orthogonal matrices are the foundation for all attitude representations.
"""

# Standard packages.
import numpy as np

# ACSToolbox packages.
from acstoolbox.foundation.math import is_column, vX


def AxisAngle(a, phi):
    """Construct a rotation matrix using the axis (a) and angle (phi) of rotation."""
    is_column(a)
    return (
        (np.cos(phi) * np.identity(3))
        + (1 - np.cos(phi)) * np.outer(a, a)
        - np.sin(phi) * vX(a)
    )


def C(a, phi_deg):
    """Construct a principal rotation. The matrix C represents
    an affine transformation from one frame to another which is
    offset by a rotation about a principal axis (a) in the first
    frame through an angle phi.

    Input:
       1. a - integer value of principal axis
       2. phi_deg - rotation angle in [degrees]
    """
    phi_rad = np.radians(phi_deg)
    c = np.cos(phi_rad)
    s = np.sin(phi_rad)

    # Construct transformation about a principal axis (1,2,3).
    C = np.zeros((3, 3))
    if a == 1:
        C = np.array([[1, 0, 0], [0, c, s], [0, -s, c]])
    elif a == 2:
        C = np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])
    elif a == 3:
        C = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    else:
        assert False, "ACSToolbox/Foundation *Error*: Principal axis out of bounds"

    return C
