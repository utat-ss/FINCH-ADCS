# Standard libraries.
import numpy as np
import math

# ACS Toolbox.
from acs_toolbox.foundation import so3 as C


class AlignConstrained:
    def __init__(self, p_b, s_b):
        """Align Constrained Trajectory
        Constructor: 1) p_b: body-fixed vector to align with
                     2) s_b: body-fixed vector to constrain with
        Outputs:     Inertial trajectory which *aligns* p_b with an inertial
                     vector a_i. Since a simple alignment results in an indeterminate
                     attitude, it is *constrained* by rotating about p_b to
                     align s_b as closely as possible with an inertial vector c_i.
        TODO:        1) Attitude: Protect against +/- infinity (+/- pi/2 rotation) in atan comp
        """

        # Normalize the align/constrain vectors.
        self.p_b_ = p_b / np.linalg.norm(p_b)
        self.s_b_ = s_b / np.linalg.norm(s_b)

    def Attitude(self, a_i, c_i):
        """
        Input: 1) a_i: Inertial vector which p_b aligns with
               2) c_i: Inertial vector which p_b constrains with
        """

        # Normalize vectors.
        a_i = a_i / np.linalg.norm(a_i)
        c_i = c_i / np.linalg.norm(c_i)

        # Rotation from inertial to aligned.
        a_a = np.cross(self.p_b_, a_i)
        phi_ai = math.acos(self.p_b_.dot(a_i) / np.linalg.norm(a_i))
        C_ai = C.AxisAngle(a_a, phi_ai)

        # Rotation from aligned to constrained.
        c_a = C_ai.dot(c_i)
        ppt_mone = np.outer(self.p_b_, self.p_b_) - np.identity(3)
        phi_ca = math.atan(
            self.s_b_.dot(np.cross(self.p_b_, c_a)) / self.s_b_.dot(ppt_mone.dot(c_a))
        )
        C_ca = C.AxisAngle(self.p_b_, phi_ca)

        return C_ai.dot(C_ca)
