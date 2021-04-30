## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
##                         ___   ___  ________    ________         ______    ________ ________                        ##
##                        |  |  |  | |   ____|   /       |        /  __  \  |   ____||   ____|                        ##
##                        |  |__|  | |  |__     |   (----` ______|  |  |  | |  |__   |  |__                           ##
##                        |   __   | |   __|     \   \    |______|  |  |  | |   __|  |   __|                          ##
##                        |  |  |  | |  |____.----)   |          |  `--'  | |  |     |  |                             ##
##                        |__|  |__| |_______|_______/            \______/  |__|     |__|                             ##
##                                                                                                                    ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##

# Import packages
import pdb
import time
import copy
import os.path
import numpy as np
import matplotlib.pyplot as plt

from scipy.io            import loadmat
from scipy.optimize      import root_scalar
from scipy.interpolate   import interp1d
from importlib_resources import files


from . import elgrid_models
from . import process_models
from . import combustion_models

# Define font settings
fontsize = 12
plt.rc('text', usetex=False)
plt.rcParams['font.family']      = 'serif'                # 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'
plt.rcParams['font.serif']       = 'times new roman'      # 'cmr10', 'palatino', 'times new roman'
plt.rcParams['mathtext.fontset'] = 'stix'                 # 'cm' (latex style), 'stix' (times new roman style), 'stixsans'



## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##




class IntegratedModel:
    """
    The IntegratedModel object defines and stores parameters of the HES-OFF concept
    """

    def __init__(self, IN):

        # Declare input variables as instance variables
        self.IN             = copy.deepcopy(IN)

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
        self.WIND_DATA = process_models.read_wind_data(IN["WIND_FILENAME"])
        self.WIND_SPEED = self.WIND_DATA["speed"]
        self.WIND_TIME = self.WIND_DATA["time"]

        # Initialize components
        self.GT = process_models.GT(model=self.GT_MODEL, number_of_units=self.GT_UNITS)
        self.WT = process_models.WT(model=self.WT_MODEL, rated_power=self.WT_RATED_POWER, wind_speed_factor=(self.WT_HUB_HEIGHT, self.WT_REF_HEIGHT))
        self.EL = process_models.EL(model=self.EL_MODEL, rated_power=self.EL_RATED_POWER, efficiency_coefficients=self.EL_EFFICIENCY)
        self.FC = process_models.FC(model=self.FC_MODEL, rated_power=self.FC_RATED_POWER, efficiency_coefficients=self.FC_EFFICIENCY)

    def __str__(self):
        class_info = '\nHES-OFF concept specifications\n'
        for key, value in self.IN.items():
            class_info += "{:<24}{}\n".format(key, value)
        return class_info


    def evaluate_process_model(self):

        #  Check the size of the power demand and heat demand arrays
        self.POWER_DEMAND, self.HEAT_DEMAND = np.atleast_1d(self.POWER_DEMAND), np.atleast_1d(self.HEAT_DEMAND)
        if self.POWER_DEMAND.size != self.HEAT_DEMAND.size:
            raise Exception("The number of elements of POWER_DEMAND and HEAT_DEMAND must be the same")

        # Initialize arrays to store the solution at each instance of the year
        p, n = self.POWER_DEMAND.size, self.WIND_SPEED.size
        flag = np.empty((p, n))  # Type of operational strategy at each time instance
        GT_power = np.empty((p, n))  # Power generated in the gas turbines (W)
        WT_power = np.empty((p, n))  # Power generated in the wind turbines (W)
        FC_power = np.empty((p, n))  # Power generated in the fuel cell system (W)
        EL_power = np.empty((p, n))  # Power consumed in the electrolyzer system (W)
        H2_level = np.empty((p, n))  # Mass of hydrogen in the storage system (kg)
        H2_produced = np.empty((p, n))  # Mass of hydrogen produced in the electrolyzer (kg)
        H2_utilized = np.empty((p, n))  # Mass of hydrogen utilized in the fuel cell + gas turbine (kg)
        CO2_produced = np.empty((p, n))  # Mass of carbon dioxide emitted to the atmosphere (kg)
        power_deficit = np.empty((p, n))  # Power demand not satisfied (W)
        # energy_surplus = np.empty((p,n))         # Extra wind energy that is dissipated (J)

        # Initialize time array (must have the same shape as the other arrays to return within dictionary)
        times = np.empty((p, n), dtype=np.int32)
        for i in range(p):
            times[i, :] = np.arange(0, n)

        # Create the natural gas + hydrogen mixture used in the gas turbines
        blend_NG_H2 = process_models.create_fluid_mixture(fluids=[combustion_models.hydrogen, combustion_models.natural_gas],
                                                   fractions=[self.GT_MAX_H2, 1 - self.GT_MAX_H2],
                                                   fraction_type="molar")

        # Loop over the time periods
        for p, (power_demand, heat_demand) in enumerate(zip(self.POWER_DEMAND, self.HEAT_DEMAND)):

            # Compute the initial level of hydrogen in the storage system
            H2_level[p, 0] = self.H2_CAPACITY * self.H2_INITIAL_LEVEL

            # Compute the minimum GT load required to satisfy the heat demand
            GT_power_min = self.GT.compute_power_from_heat(heat_demand)

            # Compute the maximum GT load of the current GT model
            GT_power_max = self.GT.rated_power

            # Compute wind power over the year
            WT_power_available = self.WT.compute_power_output(wind_speed=self.WIND_SPEED)

            for t in times[p]:

                # Use a run-out-of-steam strategy to supply the power demand
                if GT_power_min + WT_power_available[t] >= power_demand:

                    # Case 1: Use GT_min and WT to satisfy the power demand (use EL to recharge H2)
                    if H2_level[p, t] < self.H2_CAPACITY:
                        flag_current = 1
                        GT_power_current = GT_power_min
                        EL_power_current = np.minimum(self.EL_RATED_POWER, GT_power_min + WT_power_available[t] - power_demand)
                        WT_power_current = power_demand + EL_power_current - GT_power_current
                        FC_power_current = 0.00

                    # Case 2: Use GT_min and WT to satisfy the power demand (do not use EL to recharge H2)
                    else:
                        flag_current = 2
                        GT_power_current = GT_power_min
                        WT_power_current = power_demand - GT_power_current
                        EL_power_current = 0.00
                        FC_power_current = 0.00

                elif GT_power_max + WT_power_available[t] >= power_demand:

                    # Case 3: Use GT and WT to satisfy the power demand (use GT+EL to recharge H2)
                    if H2_level[p, t] < self.H2_RECHARGE_THRESHOLD * self.H2_CAPACITY:
                        flag_current = 3
                        WT_power_current = WT_power_available[t]
                        EL_power_current = np.minimum(self.EL_RATED_POWER, GT_power_max + WT_power_current - power_demand)
                        GT_power_current = np.minimum(GT_power_max, power_demand + EL_power_current - WT_power_current)
                        FC_power_current = 0.00

                    # Case 4: Use GT and WT to satisfy the power demand (do not use GT+EL to recharge H2)
                    else:
                        flag_current = 4
                        WT_power_current = WT_power_available[t]
                        GT_power_current = power_demand - WT_power_current
                        EL_power_current = 0.00
                        FC_power_current = 0.00

                else:

                    # Case 5: Use GT, WT and FC to satisfy the power demand (there is hydrogen available)
                    if H2_level[p, t] > 0.00:
                        flag_current = 5
                        WT_power_current = WT_power_available[t]
                        GT_power_current = GT_power_max
                        FC_power_current = power_demand - WT_power_current - GT_power_current
                        EL_power_current = 0.00

                    # Case 6: The GT and WT cannot satisfy the power demand (there is no hydrogen available)
                    else:
                        flag_current = 6
                        WT_power_current = WT_power_available[t]
                        GT_power_current = GT_power_max
                        FC_power_current = 0.00
                        EL_power_current = 0.00

                # Determine the type of fuel used in the gas turbines
                use_fuel_blend = H2_level[p,t] > self.H2_COFIRE_THRESHOLD * self.H2_CAPACITY
                if use_fuel_blend:
                    GT_fuel = blend_NG_H2
                else:
                    GT_fuel = combustion_models.natural_gas

                # Compute the carbon dioxide emissions (kg/s)
                CO2_produced_current = self.GT.compute_carbon_dioxide_emissions(GT_power_current, fuel=GT_fuel)*3600

                # Compute the mass flow rate of hydrogen co-fired in the gas turbine (kg/s)
                H2_cofired_current = self.GT.compute_hydrogen_mass_flow_rate(GT_power_current, fuel=GT_fuel)

                # Compute mass flow rate of hydrogen fed to the fuel cell system (kg/s)
                H2_utilized_current = self.FC.compute_hydrogen_consumption(FC_power_current)

                # Compute the mass flow rate of hydrogen produced in the electrolyzer system (kg/s)
                H2_produced_current = self.EL.compute_hydrogen_production(EL_power_current)

                # Evaluate the power balance (W)
                power_deficit_current = power_demand + EL_power_current - WT_power_current - GT_power_current - FC_power_current

                # Compute the hydrogen level for the next time instance (skip last time step computation)
                if t < times[p, -1]:
                    H2_level[p, t + 1] = H2_level[p, t] + (
                                H2_produced_current - H2_utilized_current - H2_cofired_current) * 3600

                # Store the current solution in its corresponding array
                flag[p, t] = flag_current
                GT_power[p, t] = GT_power_current
                WT_power[p, t] = WT_power_current
                FC_power[p, t] = FC_power_current
                EL_power[p, t] = EL_power_current
                power_deficit[p, t] = power_deficit_current
                CO2_produced[p, t] = CO2_produced_current
                H2_produced[p, t] = H2_produced_current
                H2_utilized[p, t] = H2_utilized_current + H2_cofired_current

        # Store the results in a dictionary
        result_dict = {"flag": flag * 1.00,  # Conversion from integer to float for Numba
                       "times": times * 1.00,  # Conversion from integer to float for Numba
                       "GT_power": GT_power,
                       "WT_power": WT_power,
                       "FC_power": FC_power,
                       "EL_power": EL_power,
                       "H2_produced": H2_produced,
                       "H2_utilized": H2_utilized,
                       "H2_level": H2_level,
                       "CO2_produced": CO2_produced,
                       "power_deficit": power_deficit}

        self.process_output = result_dict

        # Compute the accumulated CO2 emissions
        self.CO2_emissions = np.sum(self.process_output["CO2_produced"], axis=1)
        self.CO2_emissions_total = np.sum(self.CO2_emissions * self.STAGE_LENGTH)

        # Compute the accumulated energy deficit
        self.energy_deficit = np.sum(self.process_output["power_deficit"]*3600, axis=1)
        self.energy_deficit_total = np.sum(self.energy_deficit*self.STAGE_LENGTH)

    # ------------------------------------------------------------------------------------------------------------------ ##
    # Results plotting functions
    # ------------------------------------------------------------------------------------------------------------------ ##

    def plot_hydrogen_level(self):
        """ Plot hydrogen storage level over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('Hydrogen storage level over the year (kg)', fontsize=fontsize + 1, fontweight='normal', color='k')
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
            ax.legend(ncol=1, loc='lower right', fontsize=fontsize - 1, edgecolor='k', framealpha=1.0, handlelength=0.0)
            dy = np.max(y)
            ax.set_ylim([-dy / 5, np.max(y) + dy / 5])
        fig.tight_layout()
        return fig, axes

    def plot_hydrogen_balance(self):
        """ Plot the hydrogen balance over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('Hydrogen production and utilization over the year', fontsize=fontsize + 1, fontweight='normal',
                     color='k')
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
            ax.legend(ncol=1, loc='lower right', fontsize=fontsize - 1, edgecolor='k', framealpha=1.0)
            dy = max(np.max(y1) - np.min(y2), 0.02)
            ax.set_ylim([np.min(y2) - dy / 5, np.max(y1) + dy / 5])
        fig.tight_layout()
        return fig, axes

    def plot_power_deficit(self):
        """ Plot the energy deficit over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('Power deficit over the year', fontsize=fontsize + 1, fontweight='normal', color='k')
        axes = fig.subplots(n_axes)
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, self.process_output["power_deficit"][index, :] / 1e6
            for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            ax.plot(x, y, linewidth=0.75, linestyle='-', color='k', label="Period " + str(index + 1), marker="")
            ax.set_ylabel('Deficit (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
            if index + 1 == n_axes:
                ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.legend(ncol=1, loc='lower right', fontsize=fontsize - 1, edgecolor='k', framealpha=1.0, handlelength=0.0)
            dy = max(1, np.max(y))
            ax.set_ylim([-dy / 5, np.max(y) + dy / 5])
        fig.tight_layout()
        return fig, axes

    def plot_carbon_dioxide_emissions(self):
        """ Plot the carbon dioxide emissions over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('CO$_2$ emissions accumulated over the year', fontsize=fontsize + 1, fontweight='normal',
                     color='k')
        axes = fig.subplots(n_axes)
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, np.cumsum(
                self.process_output["CO2_produced"][index, :]) / 1e6
            for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            ax.plot(x, y, linewidth=0.75, linestyle='-', color='k', label="Period " + str(index + 1), marker="")
            ax.set_ylabel('CO$_2$ emissions (Mt)', fontsize=fontsize, color='k', labelpad=fontsize)
            if index + 1 == n_axes:
                ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.legend(ncol=1, loc='lower right', fontsize=fontsize - 1, edgecolor='k', framealpha=1.0, handlelength=0.0)
            dy = np.max(y)
            ax.set_ylim([-dy / 5, np.max(y) + dy / 5])
        fig.tight_layout()
        return fig, axes

    def plot_flag(self):
        """ Plot the operation flag over time """
        n_axes = self.process_output["times"].shape[0]
        fig = plt.figure(figsize=(6.0, 5.5))
        fig.suptitle('Process flag over the year', fontsize=fontsize + 1, fontweight='normal', color='k')
        axes = fig.subplots(n_axes)
        for index, ax in enumerate(axes):
            x, y = self.process_output["times"][index, :] / 24, self.process_output["flag"][index, :]
            for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
            ax.plot(x, y, linewidth=0.75, linestyle='', color='k', label="Period " + str(index + 1), marker="o",
                    markerfacecolor="w", markeredgecolor="k", markersize=3.0, markeredgewidth=0.75)
            ax.set_ylabel('Flag', fontsize=fontsize, color='k', labelpad=fontsize)
            if index + 1 == n_axes:
                ax.set_xlabel('Time (days)', fontsize=fontsize, color='k', labelpad=fontsize)
            ax.legend(ncol=1, loc='lower right', fontsize=fontsize - 1, edgecolor='k', framealpha=1.0, handlelength=0.0)
            ax.set_ylim([-0.5, 6.5])
        fig.tight_layout()
        return fig, axes
