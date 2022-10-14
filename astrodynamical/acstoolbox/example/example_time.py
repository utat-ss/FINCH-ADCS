""" This is an example of transfer time

Put the following before you run the code:
import sys
sys.path.append("FILE_PATH")

My Example:
import sys
sys.path.append("C:\\Users\\khang\\OneDrive\\Music\\Documents\\ADCS\\simulation\\ADCS\\astrodynamical")
"""

import numpy as np
from acstoolbox.time.clock import Clock


def gregorian_to_julian_date(gregorian_j2000):

    # Initialize clock and convert Gregorian date.
    clock = Clock()
    return clock.GregorianToJulianDate(gregorian_j2000)


if __name__ == "__main__":
    # example 1: J2000 (January 1, 2000 12hh00mm00ss).
    gregorian_j2000 = [2000, 1, 1, 12, 0, 0]
    JD2000 = gregorian_to_julian_date(gregorian_j2000)
    # print(JD2000)