import pdb
import copy
import numpy as np
from importlib_resources import files
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.figure import  Figure

from .utilities import *
from .process_model import *

# Define font settings
fontsize = 12
font_files = font_manager.findSystemFonts(fontpaths=[files('hes_off.core').joinpath("fonts")])
for font_file in font_files:
    font_manager.fontManager.addfont(font_file)
plt.rcParams['font.family'] = 'Arial'


class IntegratedModel:
    """
    The IntegratedModel object defines and stores parameters of the HES-OFF concept
    """

    stage_labels = ["Peak years", "Midlife years", "Tail years"]

    def __init__(self, IN):

        # Declare input variables as instance variables
        self.IN= copy.deepcopy(IN)

        # Process specifications
        self.HEAT_DEMAND = np.asarray(self.IN["HEAT_DEMAND"])*1e6
        self.POWER_DEMAND = np.asarray(self.IN["POWER_DEMAND"])*1e6
        self.STAGE_LENGTH = np.asarray(self.IN["STAGE_LENGTH"])

        # Gas turbine specifications
        self.GT_MODEL = self.IN["GT_MODEL"]
        self.GT_UNITS = self.IN["GT_UNITS"]
        self.GT_MAX_H2 = self.IN["GT_MAX_H2"]/100

        # Wind farm specifications
        self.WT_MODEL = self.IN["WT_MODEL"]
        self.WT_REF_HEIGHT = IN["WT_REF_HEIGHT"]
        self.WT_HUB_HEIGHT = IN["WT_HUB_HEIGHT"]
        self.WT_RATED_POWER = self.IN["WT_RATED_POWER"]*1e6

        # Electrolizer system specifications
        self.EL_MODEL = self.IN["EL_MODEL"]
        self.EL_RATED_POWER = self.IN["EL_RATED_POWER"]*1e6
        try:
            self.EL_EFFICIENCY = np.asarray(self.IN["EL_EFFICIENCY"])
            if self.EL_EFFICIENCY[0] is None:
                self.EL_EFFICIENCY = np.zeros((1,))
        except:
            self.EL_EFFICIENCY = np.zeros((1,))

        # Fuel cell system specifications
        self.FC_MODEL = self.IN["FC_MODEL"]
        self.FC_RATED_POWER = self.IN["FC_RATED_POWER"]*1e6
        try:
            self.FC_EFFICIENCY = np.asarray(self.IN["FC_EFFICIENCY"])
            if self.FC_EFFICIENCY[0] is None:
                self.FC_EFFICIENCY = np.zeros((1,))
        except:
            self.FC_EFFICIENCY = np.zeros((1,))

        # Hydrogen storage specifications
        self.H2_CAPACITY = self.IN["H2_CAPACITY"]
        self.H2_INITIAL_LEVEL = self.IN["H2_INITIAL_LEVEL"] / 100
        self.H2_RECHARGE_THRESHOLD = self.IN["H2_RECHARGE_THRESHOLD"]/100
        self.H2_COFIRE_THRESHOLD = self.IN["H2_COFIRE_THRESHOLD"]/100

        # Wind data specifications
        self.WIND_FILENAME = self.IN["WIND_FILENAME"]
        self.WIND_DATA = read_wind_data(IN["WIND_FILENAME"])
        self.WIND_SPEED = self.WIND_DATA["speed"]
        self.WIND_TIME = self.WIND_DATA["time"]

        # Check the input variable values
        self.check_input_variables()


    def __str__(self):
        class_info = '\nHES-OFF concept specifications\n'
        for key, value in self.IN.items():
            class_info += "{:<24}{}\n".format(key, value)
        return class_info


    def check_input_variables(self):
        """ Check that the values of the input variables are feasible """

        if not (self.HEAT_DEMAND.size == self.POWER_DEMAND.size == self.STAGE_LENGTH.size):
            raise Exception("The number of elements of POWER_DEMAND, HEAT_DEMAND and STAGE_LENGTH must be the same")

        if np.any(self.HEAT_DEMAND < 0.0):
            raise Exception("HEAT_DEMAND values must be positive")

        if np.any(self.POWER_DEMAND < 0.0):
            raise Exception("POWER_DEMAND values must be positive")

        if np.any(self.STAGE_LENGTH < 0.0):
            raise Exception("STAGE_LENGTH values must be positive")

        if self.GT_MAX_H2 < 0.0 or self.GT_MAX_H2 > 1.00:
            raise Exception("GT_MAX_H2 must be between zero and one")

        if self.GT_UNITS < 0 or not self.GT_UNITS.is_integer():
            raise Exception("GT_UNITS must be a positive integer")

        if self.WT_RATED_POWER < 0.0:
            raise Exception("WT_RATED_POWER must be positive or zero")

        if self.WT_HUB_HEIGHT < 0.0:
            raise Exception("WT_HUB_HEIGHT must be positive")

        if self.WT_REF_HEIGHT < 0.0:
            raise Exception("WT_REF_HEIGHT must be positive")

        if self.WT_HUB_HEIGHT < self.WT_REF_HEIGHT:
            raise Exception("WT_HUB_HEIGHT must be larger than WT_REF_HEIGHT")

        if self.EL_RATED_POWER < 0.0:
            raise Exception("EL_RATED_POWER must be positive or zero")

        if self.FC_RATED_POWER < 0.0:
            raise Exception("FC_RATED_POWER must be positive or zero")

        if self.H2_CAPACITY < 0.0:
            raise Exception("H2_CAPACITY must be positive or zero")

        if self.H2_INITIAL_LEVEL < 0.0 or self.H2_INITIAL_LEVEL > 1.00:
            raise Exception("H2_INITIAL_LEVEL must be between zero and one")

        if self.H2_RECHARGE_THRESHOLD < 0.0 or self.H2_RECHARGE_THRESHOLD > 1.00:
            raise Exception("H2_RECHARGE_THRESHOLD must be between zero and one")

        if self.H2_COFIRE_THRESHOLD < 0.0 or self.H2_COFIRE_THRESHOLD > 1.00:
            raise Exception("H2_COFIRE_THRESHOLD must be between zero and one")

        if self.H2_COFIRE_THRESHOLD <= self.H2_RECHARGE_THRESHOLD:
            raise Exception("H2_RECHARGE_THRESHOLD must be lower than H2_COFIRE_THRESHOLD")

        return True

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



        # output_dict = {"CO2_emissions":}
        self.GT_energy = self.create_entry("GT_power")
        self.WT_energy = self.create_entry("WT_power")
        self.FC_energy = self.create_entry("FC_power")
        self.EL_energy = self.create_entry("EL_power")
        self.H2_utilized = self.create_entry("H2_utilized")
        self.NG_utilized = self.create_entry("NG_utilized")
        self.CO2_emissions = self.create_entry("CO2_emissions")
        self.energy_deficit = self.create_entry("power_deficit")


    def create_entry(self, name):
        value = np.sum(self.process_output[name], axis=1)*self.STAGE_LENGTH  # Units (kg)
        value = np.concatenate((value, np.sum(value)[np.newaxis]))
        return value


    def get_H2_depletion_time_GT(self):
        GT_fuel = create_fluid_mixture(fluids=(natural_gas, hydrogen), fractions=(1.00 - self.GT_MAX_H2, self.GT_MAX_H2), fraction_type="molar")
        GT_power_max = compute_GT_maximum_power(model=self.GT_MODEL, number_of_units=self.GT_UNITS)
        GT_efficiency = compute_GT_efficiency(model=self.GT_MODEL, number_of_units=self.GT_UNITS, power_output=GT_power_max)
        GT_fuel_flow = compute_GT_fuel_consumption(power_output=GT_power_max, conversion_efficiency=GT_efficiency, fuel=GT_fuel)[0]
        H2_mass_flow = GT_fuel_flow * GT_fuel["y"][-1]  # Hydrogen is the last component
        depletion_time = self.H2_CAPACITY / H2_mass_flow / 3600
        return depletion_time

    def get_H2_depletion_time_FC(self):

        # Compute mass flow rate of hydrogen fed to the fuel cell system (kg/s)
        H2_mass_flow = compute_FC_hydrogen_consumption(model=self.FC_MODEL,
                                                              efficiency_coefficients=self.FC_EFFICIENCY,
                                                              rated_power=self.FC_RATED_POWER,
                                                              power_output=self.FC_RATED_POWER)[0]

        # Compute the time required to deplete the entire hydrogen storage
        depletion_time = self.H2_CAPACITY / H2_mass_flow / 3600

        return depletion_time # Units (h)

    def get_H2_refilling_time_EL(self):

        # Compute mass flow rate of hydrogen fed to the fuel cell system (kg/s)
        H2_mass_flow = compute_EL_hydrogen_production(model=self.EL_MODEL,
                                                       efficiency_coefficients=self.EL_EFFICIENCY,
                                                       rated_power=self.EL_RATED_POWER,
                                                       power_input=self.EL_RATED_POWER)[0]

        # Compute the time required to deplete the entire hydrogen storage
        refill_time = self.H2_CAPACITY / H2_mass_flow / 3600

        return refill_time  # Units (h)


# ------------------------------------------------------------------------------------------------------------------ ##
# Plot results
# ------------------------------------------------------------------------------------------------------------------ ##

    def plot_wind_timeseries(self, is_sorted=False):
        # Plot the wind speed profile as a function of time
        fig = Figure(figsize=(6.0, 4.0))
        ax = fig.subplots(1)
        ax.set_xlabel('Time (hours)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Wind speed (m/s)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        if is_sorted:
            ax.plot(self.WIND_TIME, np.sort(self.WIND_SPEED)[::-1], linewidth=0.75, linestyle='-', color='k',
                    marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w',
                    label=None)
        else:
            ax.plot(self.WIND_TIME, self.WIND_SPEED, linewidth=0.75, linestyle='-', color='k',
                    marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w',
                    label=None)
        fig.tight_layout()
        return fig


    def plot_carbon_dioxide_emissions(self):
        """ Plot the carbon dioxide emissions over time """
        fig = Figure(figsize=(6.0, 6.5))
        fig.suptitle('CO$_2$ emissions accumulated over a year (kton)', fontsize=fontsize+1, fontweight='normal', color='k')
        axes = fig.subplots(len(self.stage_labels))
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, np.cumsum(self.process_output["CO2_emissions"][index, :]) / 1e3 / 1e3
            for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            ax.plot(x, y, linewidth=0.75, linestyle='-', color='k', label="Period "+str(index+1), marker="")
            # ax.set_ylabel('CO$_2$ emissions (Mt)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.set_ylabel(self.stage_labels[index], fontsize=fontsize, color='k', labelpad=fontsize)
            if index + 1 == len(self.stage_labels):
                ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
            # ax.legend(ncol=1, loc='lower right', fontsize=fontsize-1, edgecolor='k', framealpha=1.0, handlelength=0.0)
            dy = np.max(y)
            ax.set_ylim([-dy/5, np.max(y)+dy/5])
        fig.tight_layout()
        return fig


    def plot_power_deficit(self):
        """ Plot the energy deficit over time """
        fig = Figure(figsize=(6.0, 6.5))
        fig.suptitle('Power deficit over a year (MW)', fontsize=fontsize + 1, fontweight='normal', color='k')
        axes = fig.subplots(len(self.stage_labels))
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, self.process_output["power_deficit"][index, :] / 1e6
            for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            ax.plot(x, y, linewidth=0.75, linestyle='-', color='k', label=self.stage_labels[index], marker="")
            # ax.set_ylabel('Deficit (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.set_ylabel(self.stage_labels[index], fontsize=fontsize, color='k', labelpad=fontsize)
            if index + 1 == len(self.stage_labels):
                ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
            # ax.legend(ncol=1, loc='lower right', fontsize=fontsize - 1, edgecolor='k', framealpha=1.0, handlelength=0.0)
            dy = max(1, np.max(y))
            ax.set_ylim([-dy / 5, np.max(y) + dy / 5])
        fig.tight_layout()
        return fig


    def plot_hydrogen_level(self):
        """ Plot hydrogen storage level over time """
        fig = Figure(figsize=(6.0, 6.5))
        fig.suptitle('Hydrogen storage level over a year (kg)', fontsize=fontsize+1, fontweight='normal', color='k')
        axes = fig.subplots(len(self.stage_labels))
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, self.process_output["H2_level"][index, :]
            for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            ax.plot([0.0], [0.0], linestyle="", marker="", label="Period " + str(index + 1))
            ax.plot(x, y, linewidth=0.75, linestyle='-', color='k', label="", marker="")
            # ax.set_ylabel('H$_2$ level (kg)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.set_ylabel(self.stage_labels[index], fontsize=fontsize, color='k', labelpad=fontsize)
            if index + 1 == len(self.stage_labels):
                ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
            # ax.legend(ncol=1, loc='lower right', fontsize=fontsize-1, edgecolor='k', framealpha=1.0, handlelength=0.0)
            dy = np.max(y)
            dy = np.maximum(dy, 1.00)
            ax.set_ylim([-dy/5, np.max(y)+dy/5])
        fig.tight_layout()
        return fig


    def plot_power_balance(self):
        """ Plot hydrogen storage level over time """
        fig = Figure(figsize=(6.0, 6.5))
        fig.suptitle('Power balance (MW)', fontsize=fontsize + 1, fontweight='normal', color='k')
        axes = fig.subplots(len(self.stage_labels))
        for index, ax in enumerate(axes):
            x = self.process_output["times"][index, :] / 24
            y1 = self.process_output["GT_power"][index, :] / 1e6
            y2 = self.process_output["WT_power"][index, :] / 1e6
            y3 = self.process_output["FC_power"][index, :] / 1e6
            y4 = self.process_output["EL_power"][index, :] / 1e6
            y5 = y1 + y2 + y3 - y4
            for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            ax.stackplot(x, y1, y2, y3, labels=["GT", "WT", "FC"], colors=["orangered", "forestgreen", "black"])
            ax.stackplot(x, -y4, labels=["EL"], colors=["royalblue"])
            # ax.set_ylabel('Power (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.set_ylabel(self.stage_labels[index], fontsize=fontsize, color='k', labelpad=fontsize)
            if index + 1 == len(self.stage_labels):
                ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.legend(ncol=1,bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=fontsize - 1, edgecolor='k', framealpha=1.0)
            # dy = np.max(y)
            # dy = np.maximum(dy, 1.00)
            # ax.set_ylim([-dy / 5, np.max(y) + dy / 5])
            ax.set_xlim([0, x[-1]])
        fig.tight_layout()
        return fig


    def plot_hydrogen_balance(self):
        """ Plot the hydrogen balance over time """
        fig = Figure(figsize=(6.0, 6.5))
        fig.suptitle('Hydrogen balance over a year (kg/s)', fontsize=fontsize+1, fontweight='normal', color='k')
        axes = fig.subplots(len(self.stage_labels))
        for index, ax in enumerate(axes):
            x = self.process_output["times"][index, :] / 24
            y1 = self.process_output["H2_cofired"][index, :]
            y2 = self.process_output["H2_utilized"][index, :]
            y3 = self.process_output["H2_produced"][index, :]
            for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            ax.stackplot(x, y3, labels=["EL"], colors=["royalblue"])
            ax.stackplot(x, -y1, -y2, labels=["GT", "FC"], colors=["orangered", "black"])
            # ax.set_ylabel('Mass flow (kg/s)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.set_ylabel(self.stage_labels[index], fontsize=fontsize, color='k', labelpad=fontsize)
            if index + 1 == len(self.stage_labels):
                ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.legend(ncol=1, title="Period " + str(index + 1), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=fontsize - 1, edgecolor='k', framealpha=1.0)
            dy = np.max([np.max(y1+y2), np.max(y3), 0.002])
            ax.set_ylim([-1*dy, 1*dy])
        fig.tight_layout()
        return fig


    def plot_flag(self):
        """ Plot the operation flag over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('Process flag over a year', fontsize=fontsize+1, fontweight='normal', color='k')
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
        return fig


    def plot_sensitivity_analysis(self, variable, low_value, high_value):

        # Create a copy of the input dictionary
        IN = copy.deepcopy(self.IN)

        # Compute the performance for a range of input variables
        values = np.linspace(low_value, high_value, 21)
        CO2_emissions, energy_deficit = [], []
        print("Starting", variable, "sensitivity analysis...")
        for i, value in enumerate(values):
            print_progress(i, len(values))
            IN[variable] = value
            EnergySystem = IntegratedModel(IN)
            EnergySystem.evaluate_process_model()
            CO2_emissions.append(EnergySystem.CO2_emissions_total)
            energy_deficit.append(EnergySystem.energy_deficit_total)

        # Convert list to Numpy arrays
        CO2_emissions = np.asarray(CO2_emissions)
        energy_deficit = np.asarray(energy_deficit)

        # Plot the results
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        ax.set_xlabel(variable, fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('CO$_2$ emissions (Mt)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        ax.plot(values, CO2_emissions/1e6, linewidth=0.75, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()

        # Plot the results
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        ax.set_xlabel(variable, fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Energy deficit (MWh)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        ax.plot(values, energy_deficit, linewidth=0.75, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()

        return fig


