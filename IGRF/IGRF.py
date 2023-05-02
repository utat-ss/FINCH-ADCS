import os
import numpy as np

# Class to put the igrf file values into
class igrf:
    def __init__(self, time, coeffs, parameters, terms):
        self.time = time
        self.coeffs = coeffs
        self.parameters = parameters
        self.terms = terms


filepath = "IGRF13.txt"


def readIGRF13():
    # check correct
    with open(filepath, "r") as f:

        data = np.array([])
        for line in f.readlines():

            if line[0] == "#":
                continue
            read_line = np.fromstring(line, sep=" ")
            if read_line.size == 7:
                name = os.path.split(filepath)[1]  # file name string
                values = [name] + read_line.astype(int).tolist()
            else:
                data = np.append(data, read_line)
        # unpack parameter line
        keys = ["TXT", "nmin", "nmax", "N", "order", "step", "start_year", "end_year"]
        parameters = dict(zip(keys, values))

        time = data[: parameters["N"]]
        coeffs = data[parameters["N"] :].reshape((-1, parameters["N"] + 2))
        terms = np.array(coeffs[0:, :2])
        coeffs = np.squeeze(coeffs[:, 2:])  # discard columns with n and m
    return igrf(time, coeffs, parameters, terms)


def Kronecherdelta(m, n):
    # should not need to check
    if m == n:
        return 1
    else:
        return 0


def Schmidtfactor(n, m):
    # TODO ---> double check
    if n == 0 and m == 0:
        return 1
    elif m == 0:
        return Schmidtfactor(n - 1, m) * ((2 * n - 1) / n)
    else:
        return Schmidtfactor(n, m - 1) * np.sqrt(
            ((n - m + 1) * (Kronecherdelta(m, 1) + 1)) / (n + m)
        )


def kfuction(n, m):
    # TODO ---> double check
    if n > 1:
        knm = (((n - 1) ** 2) - (m**2)) / ((2 * n - 1) * (2 * n - 3))
    elif n == 1:
        knm = 0
    return knm


def Gaussfactor(n, m, theta):
    # TODO ---> double check
    # this fuc
    if n == 0 and m == 0:
        return 1
    elif n == m:
        return np.sin(theta) * Gaussfactor(n - 1, m - 1, theta)
    elif ((n - m) % 2) != 0:
        return np.cos(theta) * Gaussfactor(n - 1, m, theta)
    else:
        return np.cos(theta) * Gaussfactor(n - 1, m, theta) - kfuction(
            n, m
        ) * Gaussfactor(n - 2, m, theta)


def Gnm(n, m, time):
    # check correct
    IGRF = readIGRF13()
    terms = IGRF.terms
    value_find = [n, m]
    x = 0
    for i in range(len(terms)):
        if np.allclose(value_find, terms[i]):
            break
        else:
            x = x + 1
    gcoeffs = IGRF.coeffs

    #print(gcoeffs[x, time])
    gnm = Schmidtfactor(n, m) * gcoeffs[x, time]
    return gnm


def hnm(n, m, time):
    # check correct
    IGRF = readIGRF13()
    terms = IGRF.terms
    value_find = [n, -m]
    x = 0
    for i in range(len(terms)):
        if np.allclose(value_find, terms[i]):
            break
        else:
            x = x + 1
    gcoeffs = IGRF.coeffs

    #print(gcoeffs[x, time])
    hnm = Schmidtfactor(n, m) * gcoeffs[x, time]
    return hnm


def gradientGaussfactor(n, m, theta):
    if n <= 0 or m == 0:
        return 0
    elif n == m:
        return np.sin(theta) * gradientGaussfactor(n - 1, m - 1, theta) + np.cos(
            theta
        ) * gradientGaussfactor(n - 1, m - 1, theta)
    elif ((n - m) % 2) != 0:
        return np.cos(theta) * gradientGaussfactor(n - 1, m, theta) - np.sin(
            theta
        ) * gradientGaussfactor(n - 1, m, theta)
    else:
        return (
            np.cos(theta) * gradientGaussfactor(n - 1, m, theta)
            - np.sin(theta) * gradientGaussfactor(n - 1, m, theta)
            - gradientGaussfactor(n - 2, m, theta) * kfuction(n, m)
        )


def yeartime(year):
    # check correct
    IGRF = readIGRF13()
    years = IGRF.time
    time = 0
    for i in range(len(years)):
        if year == years[i]:
            break
        elif year >= years[i]:
            time = time + 1
        elif year <= years[i]:

            break
    return time


def IGRFinECEF(R, theta, phi, year):
    IGRF = readIGRF13()
    time = yeartime(year)
    theta_rad = np.radians(theta)
    phi_rad = np.radians(phi)
    maxterm = IGRF.parameters["nmax"]
    a = 6371.2
    # initialize the magneticfield in each coordinate.
    B_r = 0
    B_r_in = 0
    B_theta = 0
    B_theta_in = 0
    B_phi = 0
    B_phi_in = 0
    for n in range(1,maxterm):
        if n == 0:
            break
        for m in range(n):
            B_r_in = B_r_in + (
                Gnm(n, m, time) * np.cos(m * phi_rad)
                + hnm(n, m, time) * np.sin(m * phi_rad)
            ) * Gaussfactor(n, m, theta)
            B_theta_in = B_theta_in + (
                Gnm(n, m, time) * np.cos(m * phi_rad)
                + hnm(n, m, time) * np.sin(m * phi_rad)
            ) * gradientGaussfactor(n, m, theta)
            B_phi_in = (
                m
                * (
                    (-Gnm(n, m, time)) * np.cos(m * phi_rad)
                    + hnm(n, m, time) * np.sin(m * phi_rad)
                )
                * Gaussfactor(n, m, theta)
            )
        B_r = B_r + ((a / r) ** (n + 2)) * (n + 1) * B_r_in
        B_theta = B_theta + ((a / r) ** (n + 2)) * B_theta_in
        B_phi = ((a / r) ** (n + 2)) * B_phi_in
        # reset the value
        B_r_in = 0
        b_phi_in = 0
        B_theta_in = 0
    B_phi = B_phi * ((-1) / np.sin(theta))

    return np.array([B_r, B_theta, B_phi])
