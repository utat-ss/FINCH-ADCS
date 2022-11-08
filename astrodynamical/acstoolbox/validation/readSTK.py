import os
import numpy as np

class stk:
    def __init__(self, filepath):
        self.filepath = filepath
        # check correct
        with open(self.filepath, "r") as f:
            Lines = f.readlines()
            data = []
            count = 0
            for line in Lines:
                count += 1
                no_line = line.strip()
                data.append(no_line)
            arr = []
            vector = []
            for i in data:
                arr = i.split()
                string = np.array(arr)
                floatarray = string.astype(float)
                vector.append(floatarray)

        self.vector = vector


