import h5py
import os
import matplotlib


for file in os.listdir():
    if file.startswith("out"):
        print(file)
        f = h5py.File(file, 'r')

        level0 = f.get(list(f.keys())[0])
        meas = level0.get(list(level0.keys())[0])
        magn = meas.get(list(meas.keys())[0])
        print(dict(magn.attrs.items()))
