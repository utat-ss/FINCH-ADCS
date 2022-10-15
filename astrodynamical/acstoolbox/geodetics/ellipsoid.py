""" This program will Determining Site Components"""

# acstoolbox packages.
from acstoolbox.constants.celestial_constants import *

# Standard packages.
import numpy as np


def RadiusOfCurvatureInTheMeridian(phi_d):
    N_phi = kr_E / (np.sqrt(1 - (kepsilon_e**2) * (np.sin(np.radians(phi_d)) ** 2)))
    return N_phi


def EarthSvariable(phi_d):
    S_earth = (kr_E * (1 - (kepsilon_e**2))) / (
        np.sqrt(1 - (kepsilon_e**2) * (np.sin(np.radians(phi_d)) ** 2))
    )
    return S_earth


def PossitionVectorHorizontal(phi_d, height):
    N_phi = RadiusOfCurvatureInTheMeridian(phi_d)
    Horizontal_Comp = (N_phi + height) * np.cos(np.radians(phi_d))
    return Horizontal_Comp  # Unit in KM


def PossitionVectorVertical(phi_d, height):
    S_earth = EarthSvariable(phi_d)
    Vertical_comp = (S_earth + height) * (np.sin(np.radians(phi_d)))
    return Vertical_comp  # Unit in KM


def SitePositionVector2D(phi_d, height):
    # determine site coordinate
    Vertical = PossitionVectorVertical(phi_d, height)
    Horizontal = PossitionVectorHorizontal(phi_d, height)
    # Site vector
    Site_Vector = np.array([Horizontal, Vertical])
    return Site_Vector


def SitePositionVector3D(phi_d, lamda, height):
    # determine site coordinate
    C_E = RadiusOfCurvatureInTheMeridian(phi_d)
    S_E = EarthSvariable(phi_d)

    # Calculate Site Vector
    phi_d = np.radians(phi_d)
    lamda = np.radians(lamda)
    r_i = (C_E + height) * np.cos(phi_d) * np.cos(lamda)
    r_j = (C_E + height) * np.cos(phi_d) * np.sin(lamda)
    r_k = (S_E + height) * np.sin(phi_d)
    Site_vector = np.array([r_i, r_j, r_k])
    return Site_vector
