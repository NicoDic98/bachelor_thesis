import h5py
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.special


def magnetization_exact(beta):
    if beta < 0.440686793509772:
        return 0
    else:
        return (1. - 1. / np.sinh(2 * beta) ** 4) ** (1 / 8)


def ene_exact(J):  # exact internal energy in thermodynamic limit (for h=0)\n",
    J2 = 2 * J
    k = 4 * np.sinh(J2) ** 2 / np.cosh(J2) ** 4
    eK = scipy.special.ellipk(k)
    answer = 1 + (2 / np.pi) * (2 * np.tanh(J2) ** 2 - 1) * eK
    return -J * np.cosh(J2) * answer / np.sinh(J2)


betas = []
inverse_betas = []
magnetizations = []
magnetizations_errors = []
magnetization_auto_correlation = np.zeros(10)
magnetizations_squared = []
magnetizations_squared_errors = []
energies = []
energies_errors = []
energies_squared = []
energies_squared_errors = []

for file in os.listdir():
    if file.startswith("out_02"):
        print(file)
        f = h5py.File(file, 'r')

        level0_group = f.get("level0")

        measurements_group = level0_group.get("measurements")

        betas.append(level0_group.attrs["beta"])
        inverse_betas.append(1. / betas[-1])

        if "magnetization" in measurements_group:
            magnetization_group = measurements_group.get("magnetization")
            magnetizations.append(magnetization_group.attrs["bootstrap_mean"])
            magnetizations_errors.append(np.sqrt(magnetization_group.attrs["bootstrap_variance"]))
            magnetization_auto_correlation_dataset = magnetization_group.get("auto_correlation")
            magnetization_auto_correlation = np.resize(magnetization_auto_correlation,
                                                       magnetization_auto_correlation_dataset.size)
            magnetization_auto_correlation_dataset.read_direct(magnetization_auto_correlation,
                                                               np.s_[0:magnetization_auto_correlation_dataset.size],
                                                               np.s_[0:magnetization_auto_correlation_dataset.size])

        if "magnetization_squared" in measurements_group:
            magnetization_squared_group = measurements_group.get("magnetization_squared")
            magnetizations_squared.append(magnetization_squared_group.attrs["bootstrap_mean"])
            magnetizations_squared_errors.append(np.sqrt(magnetization_squared_group.attrs["bootstrap_variance"]))

        if "energy" in measurements_group:
            energy_group = measurements_group.get("energy")
            energies.append(energy_group.attrs["bootstrap_mean"])
            energies_errors.append(np.sqrt(energy_group.attrs["bootstrap_variance"]))

        if "energy_squared" in measurements_group:
            energy_squared_group = measurements_group.get("energy_squared")
            energies_squared.append(energy_squared_group.attrs["bootstrap_mean"])
            energies_squared_errors.append(np.sqrt(energy_squared_group.attrs["bootstrap_variance"]))

fig, ax = plt.subplots()
fig: plt.Figure
ax: plt.Axes

# magnetizations
if len(magnetizations) > 0:
    ax.errorbar(inverse_betas, magnetizations, magnetizations_errors, fmt='o')
    beta_lin = np.linspace(0.25, 3, 1000)
    m_exact = np.array([magnetization_exact(temp) for temp in beta_lin])
    ax.plot(1. / beta_lin, m_exact)
    ax.plot(1. / beta_lin, -m_exact)

    ax.set_xlabel(r"1/$\beta$")
    ax.set_ylabel(r"m")
    plt.savefig("magnetisation.png")
    ax.clear()

    ax.scatter(np.arange(magnetization_auto_correlation.size),magnetization_auto_correlation)
    #ax.set_yscale("log")
    plt.savefig("magnetization_auto_correlation.png")
    ax.clear()
    ax.set_yscale("linear")

# magnetizations_squared
if len(magnetizations_squared) > 0:
    ax.errorbar(inverse_betas, magnetizations_squared, magnetizations_squared_errors, fmt='o')
    m_squared_exact = m_exact ** 2  # todo
    ax.plot(1. / beta_lin, m_squared_exact)

    ax.set_xlabel(r"1/$\beta$")
    ax.set_ylabel(r"m²")
    plt.savefig("magnetization_squared.png")
    ax.clear()

# energies
if len(energies) > 0:
    ax.errorbar(inverse_betas, energies, energies_errors, fmt='o')
    e_exact = np.array([ene_exact(temp) / temp for temp in beta_lin])
    ax.plot(1. / beta_lin, e_exact)

    ax.set_xlabel(r"1/$\beta$")
    ax.set_ylabel(r"e")
    plt.savefig("energy.png")
    ax.clear()

# energies_squared
if len(energies_squared) > 0:
    ax.errorbar(inverse_betas, energies_squared, energies_squared_errors, fmt='o')

    ax.set_xlabel(r"1/$\beta$")
    ax.set_ylabel(r"e²")
    plt.savefig("energy_squared.png")
    ax.clear()
