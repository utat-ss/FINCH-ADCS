''' This file is for making graph to validate with STK'''
import sys
sys.path.append("C:\\Users\\khang\\OneDrive\\Music\\Documents\\ADCS\\simulation")
# ACSToolbox packages.
from acstoolbox.time.clock import Clock
from acstoolbox.ephemerides.sun import Sun
from acstoolbox.validation.readSTK import stk

import matplotlib.pyplot as plt
import numpy as np
import math

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

if __name__ == "__main__":
    clock = Clock()
    sun = Sun(clock)
    # filepath = "Moon_vector_simple.txt"
    # STK = stk(filepath)


    epoch_gregorian_utc = [2018, 1, 1, 0, 0, 0]

    day = 0
    Sun_vector = []

    for i in range(1,13):
        epoch_gregorian_utc[1] = i
        if i == 1:
            for j in range(1,32):
                day = day + 1
                epoch_gregorian_utc[2] = j
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)

        if i == 2:
            for k in range(1,29):
                day = day + 1
                epoch_gregorian_utc[2] = k
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 3:
            for j in range(1,32):
                day = day + 1
                epoch_gregorian_utc[2] = j
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 4:
            for k in range(1,31):
                day = day+1
                epoch_gregorian_utc[2] = k
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 5:
            for j in range(1,32):
                day = day + 1
                epoch_gregorian_utc[2] = j
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 6:
            for k in range(1,31):
                day = day+1
                epoch_gregorian_utc[2] = k
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 7:
            for j in range(1,32):
                day = day + 1
                epoch_gregorian_utc[2] = j
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 8:
            for k in range(1,32):
                day = day+1
                epoch_gregorian_utc[2] = k
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 9:
            for j in range(1,31):
                day = day + 1
                epoch_gregorian_utc[2] = j
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 10:
            for j in range(1,32):
                day = day + 1
                epoch_gregorian_utc[2] = j
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 11:
            for k in range(1,31):
                day = day+1
                epoch_gregorian_utc[2] = k
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)
        if i == 12:
            for k in range(1,32):
                day = day+1
                epoch_gregorian_utc[2] = k
                c_tdb = clock.UTCGregotianToTUT1(epoch_gregorian_utc)
                temp = sun.GeTKmMODposittionfromUTC(c_tdb)*149597870.691
                Sun_vector.append(temp)

    print(Sun_vector)
# x_value_py = []
# for i in range(len(Sun_vector)):
#     x_value_py.append(Sun_vector[i][2])
# y_value_STK = []
# for i in range(len(Sun_vector)):
#     y_value_STK.append(STK.vector[i][2])
#
#     lst = list(range(1,366))
    # plt.plot(lst, x_value_py , 'r', label='ACStoolbox')
    # plt.plot(lst,y_value_STK, 'g',label='STK')
    # plt.title("Z value vs Day in year")
    # plt.xlabel("Day in Year")
    # plt.ylabel("Z value direction")
    # plt.legend()
    # plt.show()

    # find angle
    # sunangle = []
    # for i in range(len(Sun_vector)):
    #     anglesun =  angle(Sun_vector[i], STK.vector[i])
    #     sunangle.append(anglesun)
    # print(sunangle)
    # lst = list(range(1,366))
    # plt.plot(lst, sunangle , 'r')
    # plt.title("Angle between STK and ACStoolbox in 2018")
    # plt.xlabel("Day in Year")
    # plt.ylabel("angle in radian")
    # plt.legend()
    # plt.show()