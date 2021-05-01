import numpy as np
import matplotlib.pyplot as plt

# Define font settings
fontsize = 12
plt.rc('text', usetex=False)
plt.rcParams['font.family']      = 'serif'                # 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'
plt.rcParams['font.serif']       = 'times new roman'      # 'cmr10', 'palatino', 'times new roman'
plt.rcParams['mathtext.fontset'] = 'stix'                 # 'cm' (latex style), 'stix' (times new roman style), 'stixsans'


# ------------------------------------------------------------------------------------------------------------------ ##
# Results plotting functions
# ------------------------------------------------------------------------------------------------------------------ ##

def plot_hydrogen_level(results):
    """ Plot hydrogen storage level over time """
    n_axes = results["times"].shape[0]
    fig = plt.figure(figsize=(6.0, 5.5))
    fig.suptitle('Hydrogen storage level over the year (kg)', fontsize=fontsize+1, fontweight='normal', color='k')
    axes = fig.subplots(n_axes)
    for index, ax in enumerate(axes):
        x, y = results["times"][index, :] / 24, results["H2_level"][index, :]
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        ax.plot([0.0], [0.0], linestyle="", marker="", label="Period " + str(index + 1))
        ax.plot(x, y, linewidth=0.75, linestyle='-', color='k', label="", marker="")
        ax.set_ylabel('H$_2$ level (kg)', fontsize=fontsize, color='k', labelpad=fontsize)
        if index + 1 == n_axes:
            ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.legend(ncol=1, loc='lower right', fontsize=fontsize-1, edgecolor='k', framealpha=1.0, handlelength=0.0)
        dy = np.max(y)
        ax.set_ylim([-dy/5, np.max(y)+dy/5])
    fig.tight_layout()
    return fig, axes


def plot_hydrogen_balance(results):
    """ Plot the hydrogen balance over time """
    n_axes = results["times"].shape[0]
    fig = plt.figure(figsize=(6.0, 5.5))
    fig.suptitle('Hydrogen production and utilization over the year', fontsize=fontsize+1, fontweight='normal', color='k')
    axes = fig.subplots(n_axes)
    for index, ax in enumerate(axes):
        x1, y1 = results["times"][index, :] / 24, +results["H2_produced"][index, :]
        x2, y2 = results["times"][index, :] / 24, -results["H2_utilized"][index, :]
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        ax.plot([0.0], [0.0], linestyle="", marker="", label="Period " + str(index + 1))
        ax.plot(x1, y1, linewidth=0.75, linestyle='-', color='k', label="Produced")
        ax.plot(x2, y2, linewidth=0.75, linestyle='-', color='r', label="Utilized")
        ax.set_ylabel('Mass flow (kg/s)', fontsize=fontsize, color='k', labelpad=fontsize)
        if index + 1 == n_axes:
            ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.legend(ncol=1, loc='lower right', fontsize=fontsize-1, edgecolor='k', framealpha=1.0)
        dy = max(np.max(y1)-np.min(y2), 0.02)
        ax.set_ylim([np.min(y2)-dy/5, np.max(y1)+dy/5])
    fig.tight_layout()
    return fig, axes


def plot_power_deficit(results):
    """ Plot the energy deficit over time """
    n_axes = results["times"].shape[0]
    fig = plt.figure(figsize=(6.0, 5.5))
    fig.suptitle('Power deficit over the year', fontsize=fontsize+1, fontweight='normal', color='k')
    axes = fig.subplots(n_axes)
    for index, ax in enumerate(axes):
        x, y = results["times"][index, :] / 24, results["power_deficit"][index, :] / 1e6
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        ax.plot(x, y, linewidth=0.75, linestyle='-', color='k', label="Period "+str(index+1), marker="")
        ax.set_ylabel('Deficit (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
        if index + 1 == n_axes:
            ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.legend(ncol=1, loc='lower right', fontsize=fontsize-1, edgecolor='k', framealpha=1.0, handlelength=0.0)
        dy = max(1, np.max(y))
        ax.set_ylim([-dy/5, np.max(y)+dy/5])
    fig.tight_layout()
    return fig, axes


def plot_carbon_dioxide_emissions(results):
    """ Plot the carbon dioxide emissions over time """
    n_axes = results["times"].shape[0]
    fig = plt.figure(figsize=(6.0, 5.5))
    fig.suptitle('CO$_2$ emissions accumulated over the year', fontsize=fontsize+1, fontweight='normal', color='k')
    axes = fig.subplots(n_axes)
    for index, ax in enumerate(axes):
        x, y = results["times"][index, :] / 24, np.cumsum(results["CO2_produced"][index, :]) / 1e6
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        ax.plot(x, y, linewidth=0.75, linestyle='-', color='k', label="Period "+str(index+1), marker="")
        ax.set_ylabel('CO$_2$ emissions (Mt)', fontsize=fontsize, color='k', labelpad=fontsize)
        if index + 1 == n_axes:
            ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.legend(ncol=1, loc='lower right', fontsize=fontsize-1, edgecolor='k', framealpha=1.0, handlelength=0.0)
        dy = np.max(y)
        ax.set_ylim([-dy/5, np.max(y)+dy/5])
    fig.tight_layout()
    return fig, axes


def plot_flag(results):
    """ Plot the operation flag over time """
    n_axes = results["times"].shape[0]
    fig = plt.figure(figsize=(6.0, 5.5))
    fig.suptitle('Process flag over the year', fontsize=fontsize+1, fontweight='normal', color='k')
    axes = fig.subplots(n_axes)
    for index, ax in enumerate(axes):
        x, y = results["times"][index, :] / 24, results["flag"][index, :]
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        ax.plot(x, y, linewidth=0.75, linestyle='', color='k', label="Period "+str(index+1), marker="o",
                markerfacecolor="w", markeredgecolor="k", markersize=3.0, markeredgewidth=0.75)
        ax.set_ylabel('Flag', fontsize=fontsize, color='k', labelpad=fontsize)
        if index + 1 == n_axes:
            ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.legend(ncol=1, loc='lower right', fontsize=fontsize-1, edgecolor='k', framealpha=1.0, handlelength=0.0)
        ax.set_ylim([-0.5, 5.5])
    fig.tight_layout()
    return fig, axes
