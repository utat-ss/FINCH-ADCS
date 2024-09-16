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
def check_lat_lon_bounds(latd, latm, lond, lonm):
    # Convert to decimal degrees
    if (latd < 0):
        latm = -latm
    lat = latd + latm/60.0

    if (lond < 0):
        lonm = -lonm
    lon = lond + lonm/60.0

    return lat, lon