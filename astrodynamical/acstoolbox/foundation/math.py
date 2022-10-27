"""ACS Toolbox: Math Module
This module contains vectrix and matrix math functions.
"""

# Standard libraries.
import numpy as np


def is_column(v):
    """Check if the numpy vector is column in R3."""
    assert v.shape == (
        3,
    ), f"Vector {v} is *not* a column of size (3,). \nAssign vectors in the following format 'numpy.array([[1],[2],[3]])'"
    return 1 == 1


def vX(v):
    """Construct the cross-matrix which enables cross products in R3."""
    is_column(v)
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
