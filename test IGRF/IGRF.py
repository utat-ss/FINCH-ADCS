import os
import numpy as np

# TODO Validate this IGRF model with STK


class igrf:
    def __init__(self, filepath):
        self.filepath = filepath
        # check correct
        with open(self.filepath, "r") as f:

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
            keys = [
                "TXT",
                "nmin",
                "nmax",
                "N",
                "order",
                "step",
                "start_year",
                "end_year",
            ]
            self.parameters = dict(zip(keys, values))

            self.time = data[: self.parameters["N"]]
            coeffs = data[self.parameters["N"] :].reshape(
                (-1, self.parameters["N"] + 2)
            )
            self.terms = np.array(coeffs[0:, :2])
            self.coeffs = np.squeeze(coeffs[:, 2:])  # discard columns with n and m

    def Kronecherdelta(self, m, n):
        """ " This fuction gives the Kronecker delta value
        if n != m then Kronecker delta = 0
        if n = m then Kronecker delta = 1
        """
        if m == n:
            return 1
        else:
            return 0

    def Schmidtfactor(self, n, m):
        if n == 0 and m == 0:
            return 1
        elif m == 0:
            return self.Schmidtfactor(n - 1, m) * ((2 * n - 1) / n)
        else:
            return self.Schmidtfactor(n, m - 1) * np.sqrt(
                ((n - m + 1) * (self.Kronecherdelta(m, 1) + 1)) / (n + m)
            )

    def kfuction(self, n, m):
        if n > 1:
            knm = (((n - 1) ** 2) - (m**2)) / ((2 * n - 1) * (2 * n - 3))
        elif n == 1:
            knm = 0
        return knm

    def Gaussfactor(self, n, m, theta):
        if n == 0 and m == 0:
            return 1
        elif n == m:
            return np.sin(theta) * self.Gaussfactor(n - 1, m - 1, theta)
        elif ((n - m) % 2) != 0:
            return np.cos(theta) * self.Gaussfactor(n - 1, m, theta)
        else:
            return np.cos(theta) * self.Gaussfactor(n - 1, m, theta) - self.kfuction(
                n, m
            ) * self.Gaussfactor(n - 2, m, theta)

    def Gnm(self, n, m, time):
        terms = self.terms
        value_find = [n, m]
        x = 0
        for i in range(len(terms)):
            if np.allclose(value_find, terms[i]):
                break
            else:
                x = x + 1
        gcoeffs = self.coeffs
        gnm = self.Schmidtfactor(n, m) * gcoeffs[x, time]
        return gnm

    def hnm(self, n, m, time):
        terms = self.terms
        value_find = [n, -m]
        x = 0
        for i in range(len(terms)):
            if np.array_equal(value_find, terms[i]):
                break
            else:
                x = x + 1
        gcoeffs = self.coeffs
        hnm = self.Schmidtfactor(n, m) * gcoeffs[x, time]
        return hnm

    def gradientGaussfactor(self, n, m, theta):
        if n == 0 and m == 0:
            return 0
        elif n == m:
            return np.sin(theta) * self.gradientGaussfactor(
                n - 1, m - 1, theta
            ) + np.cos(theta) * self.gradientGaussfactor(n - 1, m - 1, theta)
        elif ((n - m) % 2) != 0:
            return np.cos(theta) * self.gradientGaussfactor(n - 1, m, theta) - np.sin(
                theta
            ) * self.gradientGaussfactor(n - 1, m, theta)
        else:
            return (
                np.cos(theta) * self.gradientGaussfactor(n - 1, m, theta)
                - np.sin(theta) * self.gradientGaussfactor(n - 1, m, theta)
                - self.gradientGaussfactor(n - 2, m, theta) * self.kfuction(n, m)
            )

    def yeartime(self, year):
        years = self.time
        time = 0
        for i in range(len(years)):
            if year == years[i]:
                break
            elif year >= years[i]:
                time = time + 1
            elif year <= years[i]:

                break
        return time

    def IGRFinECEF(self, R, theta, phi, year):

        time = self.yeartime(year)
        theta_rad = np.radians(theta)
        phi_rad = np.radians(phi)
        maxterm = self.parameters["nmax"]
        a = 6371.2
        # initialize the magneticfield in each coordinate.
        B_r = 0
        B_r_in = 0
        B_theta = 0
        B_theta_in = 0
        B_phi = 0
        B_phi_in = 0
        for n in range(1, maxterm):
            if n == 0:
                break
            for m in range(n):
                B_r_in = B_r_in + (
                    self.Gnm(n, m, time) * np.cos(m * phi_rad)
                    + self.hnm(n, m, time) * np.sin(m * phi_rad)
                ) * self.Gaussfactor(n, m, theta)
                B_theta_in = B_theta_in + (
                    self.Gnm(n, m, time) * np.cos(m * phi_rad)
                    + self.hnm(n, m, time) * np.sin(m * phi_rad)
                ) * self.gradientGaussfactor(n, m, theta)
                B_phi_in = (
                    m
                    * (
                        (-self.Gnm(n, m, time)) * np.cos(m * phi_rad)
                        + self.hnm(n, m, time) * np.sin(m * phi_rad)
                    )
                    * self.Gaussfactor(n, m, theta)
                )
            B_r = B_r + ((a / R) ** (n + 2)) * (n + 1) * B_r_in
            B_theta = B_theta + ((a / R) ** (n + 2)) * B_theta_in
            B_phi = ((a / R) ** (n + 2)) * B_phi_in
            # reset the value
            B_r_in = 0
            b_phi_in = 0
            B_theta_in = 0
        B_phi = B_phi * ((-1) / np.sin(theta))

        return np.array([B_r, B_theta, B_phi])
<<<<<<< Updated upstream
=======

    # def IGRFinECI(self, position_vector, year):
    #     z = position_vector[2]
    #     x = position_vector[0]
    #     y = position_vector[1]
    #     ab = np.linalg.norm(position_vector)
    #     # Latitude (no J2 effect)
    #     delta = np.arcsin(z/ab)
    #     #Coelevation
    #     theta = (np.pi/2) - delta
>>>>>>> Stashed changes
