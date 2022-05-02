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


def append_observable(name, measurements_group_, observable_list, observable_error_list):
    observable_group = measurements_group_.get(name)
    observable_list.append(observable_group.attrs["bootstrap_mean"])
    observable_error_list.append(np.sqrt(observable_group.attrs["bootstrap_variance"]))


def make_auto_correlation_plot_to_ax(name, measurements_group_, ax_):
    observable_group = measurements_group_.get(name)
    observable_auto_correlation_dataset = observable_group.get("auto_correlation")
    observable_auto_correlation = np.zeros(observable_auto_correlation_dataset.size)
    observable_auto_correlation_dataset.read_direct(observable_auto_correlation,
                                                    np.s_[0:observable_auto_correlation_dataset.size],
                                                    np.s_[0:observable_auto_correlation_dataset.size])

    ax_.scatter(np.arange(observable_auto_correlation.size), observable_auto_correlation)


def make_auto_correlation_plot(name, measurements_group_):
    fig_, ax_ = plt.subplots()
    fig_: plt.Figure
    ax_: plt.Axes
    make_auto_correlation_plot_to_ax(name, measurements_group_, ax_)
    ax_.set_xlabel(r"t")
    ax_.set_ylabel(r"$\Gamma$("+name+")")
    ax_.set_yscale("log")
    fig_.set_tight_layout(True)
    return fig_, ax_


def make_observable_plot(name, observable_list, observable_error_list):
    fig_, ax_ = plt.subplots()
    fig_: plt.Figure
    ax_: plt.Axes
    ax_.errorbar(inverse_betas, observable_list, observable_error_list, fmt='o')

    ax_.set_xlabel(r"1/$\beta$")
    ax_.set_ylabel(name)
    fig_.set_tight_layout(True)
    return fig_, ax_


magnetization_name = "magnetization"
magnetization_squared_name = "magnetization_squared"
energy_name = "energy"
energy_squared_name = "energy_squared"

betas = []
inverse_betas = []
magnetizations = []
magnetizations_errors = []
magnetizations_squared = []
magnetizations_squared_errors = []
energies = []
energies_errors = []
energies_squared = []
energies_squared_errors = []

# sub_folder_name = "checker_board_multi_level_hmc_2_levels/"
sub_folder_name = "a/"
# sub_folder_name = "std_hmc/"

for file in os.listdir(sub_folder_name):
    if file.startswith("out_"):
        file = sub_folder_name + file
        print(file)
        f = h5py.File(file, 'r')

        level0_group = f.get("level0")

        measurements_group = level0_group.get("measurements")

        betas.append(level0_group.attrs["beta"])
        inverse_betas.append(1. / betas[-1])

        if magnetization_name in measurements_group:
            append_observable(magnetization_name, measurements_group, magnetizations, magnetizations_errors)
            fig, ax = make_auto_correlation_plot(magnetization_name, measurements_group)
            fig.savefig(sub_folder_name + magnetization_name + "_auto_correlation.png")

        if magnetization_squared_name in measurements_group:
            append_observable(magnetization_squared_name, measurements_group, magnetizations_squared,
                              magnetizations_squared_errors)
            fig, ax = make_auto_correlation_plot(magnetization_squared_name, measurements_group)
            fig.savefig(sub_folder_name + magnetization_squared_name + "_auto_correlation.png")

        if energy_name in measurements_group:
            append_observable(energy_name, measurements_group, energies, energies_errors)
            fig, ax = make_auto_correlation_plot(energy_name, measurements_group)
            fig.savefig(sub_folder_name + energy_name + "_auto_correlation.png")

        if energy_squared_name in measurements_group:
            append_observable(energy_squared_name, measurements_group, energies_squared, energies_squared_errors)
            fig, ax = make_auto_correlation_plot(energy_squared_name, measurements_group)
            fig.savefig(sub_folder_name + energy_squared_name + "_auto_correlation.png")

# magnetizations
if len(magnetizations) > 0:
    beta_lin = np.linspace(0.25, 3, 1000)
    m_exact = np.array([magnetization_exact(temp) for temp in beta_lin])
    fig, ax = make_observable_plot(magnetization_name, magnetizations, magnetizations_errors)
    ax.plot(1. / beta_lin, m_exact)
    ax.plot(1. / beta_lin, -m_exact)
    fig.savefig(sub_folder_name + magnetization_name + ".png")

# magnetizations_squared
if len(magnetizations_squared) > 0:
    m_squared_exact = m_exact ** 2  # todo
    fig, ax = make_observable_plot(magnetization_squared_name, magnetizations_squared, magnetizations_squared_errors)

    ax.plot(1. / beta_lin, m_squared_exact)
    fig.savefig(sub_folder_name + magnetization_squared_name + ".png")

# energies
if len(energies) > 0:
    e_exact = np.array([ene_exact(temp) / temp for temp in beta_lin])
    fig, ax = make_observable_plot(energy_name, energies, energies_errors)
    ax.plot(1. / beta_lin, e_exact)
    fig.savefig(sub_folder_name + energy_name + ".png")

# energies_squared
if len(energies_squared) > 0:
    fig, ax = make_observable_plot(energy_squared_name, energies_squared, energies_squared_errors)
    # todo exact solution
    fig.savefig(sub_folder_name + energy_squared_name + ".png")
