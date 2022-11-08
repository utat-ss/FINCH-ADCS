""" This file contain calculation for the moon and sun variable """
import numpy as np


def MoonMeanAnomaliesRad(c_tdb):
    mean_anomaly_moon = np.radians((134.9 + 477198.85 * c_tdb) % 360)
    return mean_anomaly_moon


def MoonMeanArgumentofLatitudeRad(c_tdb):
    mean_argument_moon = np.radians((93.3 + 483202.03 * c_tdb) % 360)
    return mean_argument_moon


def MeanElongationfromSun(c_tdb):
    moon_elongation = np.radians((117.85 + 445267.115 * c_tdb) % 360)
    return moon_elongation


def SunMeanAnomalies(c_tdb):
    mean_anomaly_sun = np.radians((357.5 + 35999.05 * c_tdb) % 360)
    return mean_anomaly_sun


def MoonLongitude(c_tdb):
    # find moon longitude
    moon_longitude = (218.32 + 481267.8813 * c_tdb) % 360
    return moon_longitude
