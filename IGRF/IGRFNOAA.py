import os
import numpy as np
from numpy import degrees, radians
from math import pi
from scipy import interpolate


class igrf:  # A simple class to put the igrf file values into
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

    def legendre_poly(self, nmax, theta):

        costh = np.cos(radians(theta))
        sinth = np.sqrt(1 - costh**2)

        Pnm = np.zeros((nmax + 1, nmax + 2) + costh.shape)
        Pnm[0, 0] = 1  # is copied into trailing dimenions
        Pnm[1, 1] = sinth  # write theta into trailing dimenions via broadcasting

        rootn = np.sqrt(np.arange(2 * nmax**2 + 1))

        # Recursion relations after Langel "The Main Field" (1987),
        # eq. (27) and Table 2 (p. 256)
        for m in range(nmax):
            Pnm_tmp = rootn[m + m + 1] * Pnm[m, m]
            Pnm[m + 1, m] = costh * Pnm_tmp

            if m > 0:
                Pnm[m + 1, m + 1] = sinth * Pnm_tmp / rootn[m + m + 2]
            for n in np.arange(m + 2, nmax + 1):
                d = n * n - m * m
                e = n + n - 1
                Pnm[n, m] = (
                    e * costh * Pnm[n - 1, m] - rootn[d - e] * Pnm[n - 2, m]
                ) / rootn[d]
        # dP(n,m) = Pnm(m,n+1) is the derivative of P(n,m) vrt. theta
        Pnm[0, 2] = -Pnm[1, 1]
        Pnm[1, 2] = Pnm[1, 0]
        for n in range(2, nmax + 1):
            Pnm[0, n + 1] = -np.sqrt((n * n + n) / 2) * Pnm[n, 1]
            Pnm[1, n + 1] = (
                np.sqrt(2 * (n * n + n)) * Pnm[n, 0]
                - np.sqrt((n * n + n - 2)) * Pnm[n, 2]
            ) / 2

            for m in np.arange(2, n):
                Pnm[m, n + 1] = 0.5 * (
                    np.sqrt((n + m) * (n - m + 1)) * Pnm[n, m - 1]
                    - np.sqrt((n + m + 1) * (n - m)) * Pnm[n, m + 1]
                )
            Pnm[n, n + 1] = np.sqrt(2 * n) * Pnm[n, n - 1] / 2
        return Pnm

    def synth_values(self, coeffs, radius, theta, phi):

        # ensure ndarray inputs
        coeffs = np.array(coeffs, dtype=float)
        radius = np.array(radius, dtype=float) / 6371.2  # Earth's average radius
        theta = np.array(theta, dtype=float)
        phi = np.array(phi, dtype=float)
        nmin = 1
        nmax = 13

        theta = theta[..., None]  # first dimension is theta
        phi = phi[None, ...]  # second dimension is phi

        # get shape of broadcasted result

        b = np.broadcast(radius, theta, phi, np.broadcast_to(0, coeffs.shape[:-1]))

        grid_shape = b.shape

        # initialize radial dependence given the source
        r_n = radius ** (-(nmin + 2))

        # compute associated Legendre polynomials as (n, m, theta-points)-array
        Pnm = self.legendre_poly(nmax, theta)

        # save sinth for fast access
        sinth = Pnm[1, 1]

        # calculate cos(m*phi) and sin(m*phi) as (m, phi-points)-array
        phi = radians(phi)
        cmp = np.cos(np.multiply.outer(np.arange(nmax + 1), phi))
        smp = np.sin(np.multiply.outer(np.arange(nmax + 1), phi))

        # allocate arrays in memory
        B_radius = np.zeros(grid_shape)
        B_theta = np.zeros(grid_shape)
        B_phi = np.zeros(grid_shape)

        num = nmin**2 - 1
        for n in range(nmin, nmax + 1):
            B_radius += (n + 1) * Pnm[n, 0] * r_n * coeffs[..., num]

            B_theta += -Pnm[0, n + 1] * r_n * coeffs[..., num]

            num += 1

            for m in range(1, n + 1):
                B_radius += (
                    (n + 1)
                    * Pnm[n, m]
                    * r_n
                    * (coeffs[..., num] * cmp[m] + coeffs[..., num + 1] * smp[m])
                )

                B_theta += (
                    -Pnm[m, n + 1]
                    * r_n
                    * (coeffs[..., num] * cmp[m] + coeffs[..., num + 1] * smp[m])
                )

                with np.errstate(divide="ignore", invalid="ignore"):
                    # handle poles using L'Hopital's rule
                    div_Pnm = np.where(theta == 0.0, Pnm[m, n + 1], Pnm[n, m] / sinth)
                    div_Pnm = np.where(theta == degrees(pi), -Pnm[m, n + 1], div_Pnm)
                B_phi += (
                    m
                    * div_Pnm
                    * r_n
                    * (coeffs[..., num] * smp[m] - coeffs[..., num + 1] * cmp[m])
                )

                num += 2
            r_n = r_n / radius  # equivalent to r_n = radius**(-(n+2))
        return B_radius, B_theta, B_phi

    def IGRFtest(self, alt, lat, lon, year):
        colat = 90 - lat
        f = interpolate.interp1d(self.time, self.coeffs, fill_value="extrapolate")
        coeffs = f(year)
        Br, Bt, Bp = self.synth_values(coeffs.T, alt, colat, lon)
        # For the SV, find the 5 year period in which the date lies and compute
        # the SV within that period. IGRF has constant SV between each 5 year period
        # We don't need to subtract 1900 but it makes it clearer:
        epoch = (year - 1900) // 5
        epoch_start = epoch * 5
        # Add 1900 back on plus 1 year to account for SV in nT per year (nT/yr):
        colat = 90 - lat
        coeffs_sv = f(1900 + epoch_start + 1) - f(1900 + epoch_start)
        Brs, Bts, Bps = self.synth_values(coeffs_sv, alt, colat, lon)

        # Use the main field coefficients from the start of each five epoch
        # to compute the SV for Dec, Inc, Hor and Total Field (F)
        # [Note: these are non-linear components of X, Y and Z so treat separately]
        coeffsm = f(1900 + epoch_start)
        Brm, Btm, Bpm = self.synth_values(coeffsm.T, alt, colat, lon)
        # Rearrange to X, Y, Z components
        X = -Bt
        Y = Bp
        Z = -Br
        # For the SV
        dX = -Bts
        dY = Bps
        dZ = -Brs
        Xm = -Btm
        Ym = Bpm
        Zm = -Brm
        magneticfield = np.array([X, Y, Z])
        return magneticfield


if __name__ == "__main__":
    R = np.linalg.norm([859.52248524, -4541.5945528, 4380.34474178])
    theta = -79.2832
    phi = 43.6532
    year = 2022
    filepath = "IGRF13.txt"
    year = 2022
    IGRF = igrf(filepath)
    # for i in range(1,5):
    #     for j in range(i):
    #         print(IGRF.Gaussfactor(i, j, theta))
    magneticfield = IGRF.IGRFtest(R,  phi, theta, year)
    print(magneticfield)