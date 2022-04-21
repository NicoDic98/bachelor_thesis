import h5py
import os
import matplotlib.pyplot as plt

betas = []
inverse_betas=[]
magnetizations = []
magnetizations_squared = []
for file in os.listdir():
    if file.startswith("out"):
        print(file)
        f = h5py.File(file, 'r')

        level0_group = f.get("level0")

        measurements_group = level0_group.get("measurements")

        betas.append(level0_group.attrs["beta"])
        inverse_betas=1./betas[-1]
        magnetization_group = measurements_group.get("magnetization")
        magnetizations.append(magnetization_group.attrs["mean"])

        magnetization_squared_group = measurements_group.get("magnetization_squared")
        magnetizations_squared.append(magnetization_squared_group.attrs["mean"])


plt.plot(betas, magnetizations)

