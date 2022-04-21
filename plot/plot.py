import h5py
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

betas = []
inverse_betas = []
magnetizations = []
magnetizations_errors = []
magnetizations_squared = []
magnetizations_squared_errors = []
for file in os.listdir():
    if file.startswith("out"):
        print(file)
        f = h5py.File(file, 'r')

        level0_group = f.get("level0")

        measurements_group = level0_group.get("measurements")

        betas.append(level0_group.attrs["beta"])
        inverse_betas.append(1. / betas[-1])
        magnetization_group = measurements_group.get("magnetization")
        magnetizations.append(magnetization_group.attrs["bootstrap_mean"])
        magnetizations_errors.append(np.sqrt(magnetization_group.attrs["bootstrap_variance"]))

        magnetization_squared_group = measurements_group.get("magnetization_squared")
        magnetizations_squared.append(magnetization_squared_group.attrs["mean"])


def magnetization_exact(beta):
    if beta < 0.440686793509772:
        return 0
    else:
        return (1. - 1. / np.sinh(2 * beta) ** 4) ** (1 / 8)


print(magnetizations_errors)
fig, ax = plt.subplots()
fig: plt.Figure
ax: plt.Axes
ax.errorbar(inverse_betas, magnetizations, magnetizations_errors, fmt='o')
beta_lin = np.linspace(0.25, 3, 1000)
m_exact = np.array([magnetization_exact(temp) for temp in beta_lin])
plt.plot(1. / beta_lin, m_exact)
plt.plot(1. / beta_lin, -m_exact)

ax.set_xlabel(r"1/$\beta$")
ax.set_ylabel(r"m")
plt.savefig("magnetisation.png")
