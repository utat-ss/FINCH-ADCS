import os
from acstoolbox.constants.time_constants import *
import math
import numpy as np
import os
import pandas as pd

# TODO
# 1. Use 61 instead of 60 in JD when the day includes a leap second (confirm for all frames)
# 1. Add check to JD when the formula becomes invalid


def ConstructEarthObservationParameterTables():
    # TODO: Confirm infer_nrows to determine negative values of columns
    # Long Term: EOP 14 C04 (IAU2000A)
    # Reference: https://hpiers.obspm.fr/eoppc/eop/eopc04/C04.guide.pdf
    eop_file = os.path.dirname(os.path.abspath(__file__)) + "/EOP.txt"
    header_list = [
        "year",
        "month",
        "date",
        "mjd",
        "x_arcsec",
        "y_arcsec",
        "dUT1_s",
        "LOD_s",
        "dX_arcsec",
        "dY_arcsec",
        "x_Err_arcsec",
        "y_err_arcsec",
        "dUT1_err_s",
        "LOD_err_s",
        "dX_err_arcsec",
        "dY_err_arcsec",
    ]
    df = pd.read_fwf(eop_file, sep=" ", names=header_list, infer_nrows=1000)

    return df


def IsLeapYear(yyyy):
    # TODO: When the year is a century-multiple,
    # it is only a leap year if divisible by 400.
    if yyyy % 4 == 0:
        return 1 == 1

    return 0 == 1


def GetDaysInEachMonth():
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    return days_in_month


def DayofYear(yyyy, mm, dd, h, m, s):
    # Day of the year is not the same as days elapsed.
    # Ex. the day of the year on Jan 1 at noon is 1.5,
    # even though a full day has not passed
    days_in_month = GetDaysInEachMonth()

    if IsLeapYear(yyyy):
        days_in_month[1] += 1

    # Numebr of days up to the month.
    # TODO: Use constants.
    days = sum(days_in_month[0:mm]) + dd + h / 24 + m / 1440 + s / 86400

    return days


def YearFractionalDaystoJDUTC(yyyy, total_fractional_days):
    days_per_month = GetDaysInEachMonth()
    if IsLeapYear(yyyy):
        days_per_month[1] = 29

    mm = 0
    full_months_in_days = 0
    for month in np.linspace(0, 11, 12):
        days_up_to_month_end = sum(days_per_month[0 : int(month + 1)])
        if days_up_to_month_end > total_fractional_days:
            mm = month + 1
            if month > 0:
                full_months_in_days = sum(days_per_month[0 : (int(month))])
            break
    print(f"full: {full_months_in_days}")
    fractional_days_remaining = total_fractional_days - full_months_in_days
    dd = math.floor(fractional_days_remaining)

    print(f"month: {mm}, day: {dd}, fract: {fractional_days_remaining-dd}")
    return UTCFractionalDayToJDUTC(yyyy, mm, dd, fractional_days_remaining - dd)


def UTCFractionalDayToJDUTC(y, m, d, fractional_day):
    return (
        367 * y
        - int(7 * (y + int((m + 9) / 12)) / 4)
        + int(275 * m / 9)
        + d
        + 1721013.5
        + fractional_day
    )


def GetdATTable():
    return {
        "51179.0": 32,
        "53736.0": 33,
        "54832.0": 34,
        "56109.0": 35,
        "57204.0": 36,
        "57754.0": 37,
    }


# TODO: Update dAT table
# TODO: Julian date for leap seconds
# TODO: Julian date for UTC vs. TT
# TODO: CLS vs SELF
# NOTE: Gregorian dates are stored as lists in the format [y m d hh mm ss]
class Clock:
    def __init__(self):
        self.eop_df_ = ConstructEarthObservationParameterTables()
        self.dAT_ = GetdATTable()

    # Evaluate the Julian Date from a Gregorian date input.
    # The Julian Date frame is the same frame as the input.
    def GregorianToJulianDate(self, gregorian):
        return (
            367 * gregorian[0]
            - int(7 * (gregorian[0] + int((gregorian[1] + 9) / 12)) / 4)
            + int(275 * gregorian[1] / 9)
            + gregorian[2]
            + 1721013.5
            + ((gregorian[5] / 60 + gregorian[4]) / 60 + gregorian[3]) / 24
        )

    # Evaluate the Julian seconds from J2000 (TT).
    def GregorianToJSJ2000(self, gregorian):
        jd_from_j2000 = self.GregorianToJulianDate(gregorian) - JD_J2000

        # This function is most useful when used in TT.
        # Recall that 1 day all time standards have 86400 s/day.
        # However, only 1 TT day and 1 TAI day contain 86400 SI seconds.
        # Other timescales (e.g., UT) employ "seconds" which depend with the Earth's variable rotation rate.
        return jd_from_j2000 * DAY_IN_SECONDS

    # Evaluate the Julian Date from Julian Seconds from J2000 (TT).
    def JSJ2000ToJD(self, js_j2000):
        # This function is most useful when used in TT.
        # Assume that each day has 86400s (see GregorianToJSJ2000)
        return (js_j2000 / DAY_IN_SECONDS) + JD_J2000

    # Evaluate Modified Julian Date.
    # The Modified Julian Date is the same frame as the JD input.
    def MJD(cls, JD):
        return JD - MODIFIED_JULIAN_DATE_OFFSET

    def JDToT(cls, jd):
        return (jd - JD_J2000) / kCenturyInJulianDays

    def UTCGregoriantoTAIJD(self, gregorian_utc):
        dat_s = self.GetdATfromGregorian(gregorian_utc)

        # Add the leap second offset to the UTC seconds.
        # Valid addition since the JD conversion divides the seconds by 60 (61) to sum minutes.
        gregorian_tai = gregorian_utc
        gregorian_tai[5] = gregorian_tai[5] + dat_s

        return self.GregorianToJulianDate(gregorian_tai)

    def TAIstoTTs(self, tai_s):
        return tai_s + ATOMIC_TO_TERRESTRIAL_S

    def UTCGregotianToTUT1(self, gregorian_utc):
        jd_utc = self.GregorianToJulianDate(gregorian_utc)
        mjd_utc = self.MJD(jd_utc)
        dut1_s = self.GetEarthObservationParameter(mjd_utc, "dUT1_s")

        # Add the UT1 offset to the UTC seconds.
        # Valid addition since the JD conversion divides the seconds by 60 (61) to sum minutes.
        gregorian_ut1 = gregorian_utc
        gregorian_ut1[5] = gregorian_ut1[5] + dut1_s
        jd_ut1 = self.GregorianToJulianDate(gregorian_ut1)

        return self.JDToT(jd_ut1)

    def GetdUT1fromGregorian(self, gregorian):
        MJD = self.MJD(self.GregorianToJulianDate(gregorian))
        return self.GetEarthObservationParameter(MJD, "dUT1_s")

    def GetEarthObservationParameter(self, mjd, param):
        mjd_lb = math.floor(mjd)
        mjd_ub = math.ceil(mjd)

        lb = pd.to_numeric(self.eop_df_.loc[self.eop_df_["mjd"] == mjd_lb, param])
        ub = pd.to_numeric(self.eop_df_.loc[self.eop_df_["mjd"] == mjd_ub, param])

        return np.interp(mjd, [mjd_lb, mjd_ub], [lb.iloc[0], ub.iloc[0]])

    def GetdATfromGregorian(self, gregorian):
        MJD = self.MJD(self.GregorianToJulianDate(gregorian))
        index = [x for x in self.dAT_.keys() if float(x) <= MJD]

        return self.dAT_[index[-1]]
