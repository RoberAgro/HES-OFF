import pdb
import copy
import numpy as np
import numba as nb
import matplotlib.pyplot as plt

from .process_model import *
from .process_model_extra import *


# Define font settings
fontsize = 12
plt.rc('text', usetex=False)
plt.rcParams['font.family']      = 'serif'                # 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'
plt.rcParams['font.serif']       = 'times new roman'      # 'cmr10', 'palatino', 'times new roman'
plt.rcParams['mathtext.fontset'] = 'stix'                 # 'cm' (latex style), 'stix' (times new roman style), 'stixsans'


class IntegratedModel:
    """
    The IntegratedModel object defines and stores parameters of the HES-OFF concept
    """

    def __init__(self, IN):

        # # Check the number of dimensions
        # if int(IN['NDIM'][0]) != 3: raise Exception("The number of dimensions must be 3!")

        # Declare input variables as instance variables
        self.IN= copy.deepcopy(IN)

        # Process specifications
        self.HEAT_DEMAND = np.asarray(self.IN["HEAT_DEMAND"])
        self.POWER_DEMAND = np.asarray(self.IN["POWER_DEMAND"])
        self.STAGE_LENGTH = np.asarray(self.IN["STAGE_LENGTH"])

        # Gas turbine specifications
        self.GT_MODEL = self.IN["GT_MODEL"]
        self.GT_UNITS = self.IN["GT_UNITS"]
        self.GT_MAX_H2 = self.IN["GT_MAX_H2"]

        # Wind farm specifications
        self.WT_MODEL = self.IN["WT_MODEL"]
        self.WT_REF_HEIGHT = IN["WT_REF_HEIGHT"]
        self.WT_HUB_HEIGHT = IN["WT_HUB_HEIGHT"]
        self.WT_RATED_POWER = self.IN["WT_RATED_POWER"]

        # Electrolizer system specifications
        self.EL_MODEL = self.IN["EL_MODEL"]
        self.EL_RATED_POWER = self.IN["EL_RATED_POWER"]
        self.EL_EFFICIENCY = np.asarray(self.IN["EL_EFFICIENCY"])

        # Fuel cell system specifications
        self.FC_MODEL = self.IN["FC_MODEL"]
        self.FC_RATED_POWER = self.IN["FC_RATED_POWER"]
        self.FC_EFFICIENCY = np.asarray(self.IN["FC_EFFICIENCY"])

        # Hydrogen storage specifications
        self.H2_CAPACITY = self.IN["H2_CAPACITY"]
        self.H2_INITIAL_LEVEL = self.IN["H2_INITIAL_LEVEL"]
        self.H2_RECHARGE_THRESHOLD = self.IN["H2_RECHARGE_THRESHOLD"]
        self.H2_COFIRE_THRESHOLD = self.IN["H2_COFIRE_THRESHOLD"]

        # Wind data specifications
        self.WIND_FILENAME = self.IN["WIND_FILENAME"]
        self.WIND_DATA = read_wind_data(IN["WIND_FILENAME"])
        self.WIND_SPEED = self.WIND_DATA["speed"]
        self.WIND_TIME = self.WIND_DATA["time"]
        
        # Results
        self.process_output = None


    def __str__(self):
        class_info = '\nHES-OFF concept specifications\n'
        for key, value in self.IN.items():
            class_info += "{:<24}{}\n".format(key, value)
        return class_info


    def evaluate_process_model(self):

        # Evaluate the process model
        self.process_output = evaluate_process_model(self.HEAT_DEMAND, self.POWER_DEMAND,
                                                     self.GT_MODEL, self.GT_UNITS, self.GT_MAX_H2,
                                                     self.WT_MODEL, self.WT_RATED_POWER,
                                                     self.WT_REF_HEIGHT, self.WT_HUB_HEIGHT,
                                                     self.EL_MODEL, self.EL_RATED_POWER, self.EL_EFFICIENCY,
                                                     self.FC_MODEL, self.FC_RATED_POWER, self.FC_EFFICIENCY,
                                                     self.H2_CAPACITY, self.H2_INITIAL_LEVEL,
                                                     self.H2_RECHARGE_THRESHOLD, self.H2_COFIRE_THRESHOLD,
                                                     self.WIND_SPEED, self.WIND_TIME,
                                                     natural_gas, hydrogen)

        # Compute the accumulated CO2 emissions
        self.CO2_emissions = np.sum(self.process_output["CO2_produced"], axis=1)
        self.CO2_emissions_total = np.sum(self.CO2_emissions * self.STAGE_LENGTH)

        # Compute the accumulated energy deficit
        self.energy_deficit = np.sum(self.process_output["power_deficit"]*3600, axis=1)
        self.energy_deficit_total = np.sum(self.energy_deficit*self.STAGE_LENGTH)




# ------------------------------------------------------------------------------------------------------------------ ##
# Wind data plotting
# ------------------------------------------------------------------------------------------------------------------ ##
    
    def plot_wind_timeseries(self, is_sorted=False, fig=None, ax=None):
        # Plot the wind speed profile as a function of time
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Time (hours)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Wind speed (m/s)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        if is_sorted:
            ax.plot(self.WIND_TIME, np.sort(self.WIND_SPEED)[::-1], linewidth=1.25, linestyle='-', color='k',
                    marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        else:
            ax.plot(self.WIND_TIME, self.WIND_SPEED, linewidth=1.25, linestyle='-', color='k',
                    marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax
    

# ------------------------------------------------------------------------------------------------------------------ ##
# Results plotting functions
# ------------------------------------------------------------------------------------------------------------------ ##
    
    def plot_hydrogen_level(self):
        """ Plot hydrogen storage level over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('Hydrogen storage level over the year (kg)', fontsize=fontsize+1, fontweight='normal', color='k')
        axes = fig.subplots(n_axes)
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, self.process_output["H2_level"][index, :]
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
    
    
    def plot_hydrogen_balance(self):
        """ Plot the hydrogen balance over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('Hydrogen production and utilization over the year', fontsize=fontsize+1, fontweight='normal', color='k')
        axes = fig.subplots(n_axes)
        for index, ax in enumerate(axes):
            x1, y1 = self.process_output["times"][index, :] / 24, +self.process_output["H2_produced"][index, :]
            x2, y2 = self.process_output["times"][index, :] / 24, -self.process_output["H2_utilized"][index, :]
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
    
    
    def plot_power_deficit(self):
        """ Plot the energy deficit over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('Power deficit over the year', fontsize=fontsize+1, fontweight='normal', color='k')
        axes = fig.subplots(n_axes)
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, self.process_output["power_deficit"][index, :] / 1e6
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
    
    
    def plot_carbon_dioxide_emissions(self):
        """ Plot the carbon dioxide emissions over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('CO$_2$ emissions accumulated over the year', fontsize=fontsize+1, fontweight='normal', color='k')
        axes = fig.subplots(n_axes)
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, np.cumsum(self.process_output["CO2_produced"][index, :]) / 1e6
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
    
    
    def plot_flag(self):
        """ Plot the operation flag over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('Process flag over the year', fontsize=fontsize+1, fontweight='normal', color='k')
        axes = fig.subplots(n_axes)
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, self.process_output["flag"][index, :]
            for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            ax.plot(x, y, linewidth=0.75, linestyle='', color='k', label="Period "+str(index+1), marker="o",
                    markerfacecolor="w", markeredgecolor="k", markersize=3.0, markeredgewidth=0.75)
            ax.set_ylabel('Flag', fontsize=fontsize, color='k', labelpad=fontsize)
            if index + 1 == n_axes:
                ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.legend(ncol=1, loc='lower right', fontsize=fontsize-1, edgecolor='k', framealpha=1.0, handlelength=0.0)
            ax.set_ylim([-0.5, 6.5])
        fig.tight_layout()
        return fig, axes
