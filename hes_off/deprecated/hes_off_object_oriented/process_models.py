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
import os.path
import numpy as np
import matplotlib.pyplot as plt

from scipy.io            import loadmat
from scipy.optimize      import root_scalar
from scipy.interpolate   import interp1d
from importlib_resources import files

from .combustion_models  import *


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


def read_wind_data(filename):

    # Load wind data file
    if filename == "SLEIPNERWIND":
        mat_data = loadmat(files('hes_off.data_files').joinpath("sleipnerwind.mat"))
    else:
        mat_data = loadmat(filename)

    # Store wind data into a dictionary
    wind_data = {"speed": mat_data["wind"][0][0][0].squeeze(),
                 "time":  mat_data["wind"][0][0][1].squeeze(),
                 "year":  mat_data["wind"][0][0][2].squeeze()}

    # Convert from days to hours
    wind_data["time"] = (wind_data["time"] - wind_data["time"][0]) * 24

    # Convert from minute-based to hourly data
    N = 60
    wind_data["time"] = wind_data["time"][0:-1:N]
    wind_data["speed"] = np.asarray([np.mean(wind_data["speed"][i:i + N]) for i in range(0, len(wind_data["speed"]), N)])

    return wind_data

def plot_wind_timeseries(wind_data, is_sorted=False, fig=None, ax=None):
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
        ax.plot(wind_data["time"], np.sort(wind_data["speed"])[::-1], linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
    else:
        ax.plot(wind_data["time"], wind_data["speed"], linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
    fig.tight_layout()
    return fig, ax


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


class WT:
    """
    The WT object defines and stores parameters of a wind turbine
    """

    def __init__(self, model=None, rated_power=None, wind_speed_factor=None):

        # Initialize gas turbine model
        if model is None:
            print("The wind turbine model was not specified. Using model='HYWIND'")
            self.model = "HYWIND"
        else:
            self.model = model

        # Initialize wind farm rated power
        if rated_power is None:
            print("The wind farm rated power was not specified. Using one single wind turbine")
            self.rated_power = self._unit_power
        else:
            self.rated_power = rated_power

        # Initialize wind speed factor
        if wind_speed_factor is None:
            print("The wind speed factor was not specified. Using wind_speed_factor=1")
            self.wind_speed_factor = 1.0
        else:
            self.wind_speed_factor = wind_speed_factor


    def __str__(self):
        class_info = '\nWind farm specifications\n'
        class_info += '{:<30}{}\n'.format("WT model:",                 self.model)
        class_info += '{:<30}{}\n'.format("WT number of units:",       self.number_of_units)
        class_info += '{:<30}{:.3f}\n'.format("WT unit power (MW):",   self.unit_power/1e6)
        class_info += '{:<30}{:.3f}\n'.format("WT rated power (MW):",  self.rated_power/1e6)
        class_info += '{:<30}{:.3f}\n'.format("WT hub height (m):",    self._hub_height)
        class_info += '{:<30}{:.3f}\n'.format("WT ref height (m):",    self._ref_height)
        class_info += '{:<30}{:.3f}\n'.format("WT wind speed factor:", self.wind_speed_factor)
        return class_info


    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model_name):
        self._model = model_name

        # Specify the power curve of each turbine model
        if model_name == 'HYWIND':
            # Based on based on https://www.uib.no/sites/w3.uib.no/files/attachments/hywind_energy_lab.pdf (page 25)
            self._unit_power = 6e6
            self._speed_data = np.asarray([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 50])
            self._power_data = np.asarray([0.0000, 0.0000, 0.0000, 0.0000, 0.0290, 0.0725, 0.1304, 0.2101, 0.3261, 0.4638,
                                           0.6232, 0.7754, 0.8913, 0.9565, 0.9855, 1.0000, 1.0000, 1.0000, 0.0000, 0.0000])
        elif model_name == 'NREL':
            # Based on report NREL/TP-500-38060
            self._unit_power = 5e6
            self._speed_data = np.asarray([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 25, 30, 100])
            self._power_data = np.asarray([0, 0, 0, 170, 391, 731, 1173, 1752, 2534, 3452, 4558, 5000, 5000, 5000,
                                           5000, 5000, 5000, 5000, 0, 0])
        else:
            raise Exception("Invalid wind turbine model\nValid options: 'HYWIND', 'NREL'")

        # Normalize power data
        self._normalized_power_data = self._power_data / np.max(self._power_data)
        self._power_data = self._normalized_power_data * self._unit_power

        # Create linear interpolant from power/speed data
        self._power_interpolant = interp1d(self._speed_data, self._normalized_power_data, kind='linear')


    @property
    def unit_power(self):
        return self._unit_power

    @property
    def number_of_units(self):
        return self.rated_power / self._unit_power

    @property
    def rated_power(self):
        return self._rated_power

    @rated_power.setter
    def rated_power(self, rated_power):
        if rated_power >= 0:
            self._rated_power = rated_power
        else:
            raise ValueError("Rated power must be positive")

    @property
    def wind_speed_factor(self):
        return self._wind_speed_factor

    @wind_speed_factor.setter
    def wind_speed_factor(self, arg):

        if np.isscalar(arg):
            self._hub_height, self._ref_height = 1.0, 1.0
            self._wind_speed_factor = arg
        elif np.size(arg) == 2:
            self._hub_height, self._ref_height = arg
            self._wind_speed_factor = (self._hub_height / self._ref_height) ** (1 / 7)
        else:
            raise Exception('Invalid wind_speed_factor specification\n'
                            'Provide a scalar value or a (hub_height, ref_height) pair')


    def compute_power_output(self, wind_speed):
        return self._power_interpolant(np.asarray(wind_speed) * self.wind_speed_factor) * self.rated_power

    def compute_energy_output(self, wind_speed, time_series):
        return np.trapz(self.compute_power_output(wind_speed), time_series * 3600)

    def compute_capacity_factor(self, wind_speed, time_series):
        return self.compute_energy_output(wind_speed, time_series) / (self.rated_power * 3600 * 24 * 365)


    def plot_power_timeseries(self, wind_speed, time_series, is_sorted=False, fig=None, ax=None):
        # Plot the wind power as a function of time
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Time (hours)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Power (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        if is_sorted:
            ax.plot(time_series, np.sort(self.compute_power_output(wind_speed))[::-1] / 1e6, linewidth=1.25, linestyle='-', color='k',
                    marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        else:
            ax.plot(time_series, self.compute_power_output(wind_speed) / 1e6, linewidth=1.25, linestyle='-', color='k',
                    marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax

    def plot_normalized_power_curve(self, fig=None, ax=None):
        # Plot the normalized power curve of the wind turbine
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Wind speed (m/s)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Normalized power (-)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        speed = np.linspace(self._speed_data[0], self._speed_data[-1], 250)
        power = self._power_interpolant(speed)
        ax.plot(speed, power, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        ax.plot(self._speed_data, self._normalized_power_data, linewidth=1.25, linestyle=' ', color='k',
                marker='o', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax

    def plot_power_curve(self, fig=None, ax=None):
        # Plot the power curve of the wind turbine
        fontsize = 14
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wind speed (m/s)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Power output (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        speed = np.linspace(self._speed_data[0], self._speed_data[-1], 250)
        power = self._power_interpolant(speed) * self.rated_power
        ax.plot(speed, power/1e6, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        ax.plot(self._speed_data, self._normalized_power_data * self.rated_power / 1e6, linewidth=1.25, linestyle=' ', color='k',
                marker='o', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


class GT:
    """
    The GT object defines and stores parameters of a gas turbine
    """

    def __init__(self, model=None, number_of_units=None):

        # Initialize gas turbine model
        if model is None:
            print("The gas turbine model was not specified. Using model='LM2500+G4'")
            self.model = "LM2500+G4"
        else:
            self.model = model

        # Initialize the number of gas turbine units
        if number_of_units is None:
            print("The number of gas turbines was not specified. Using one single gas turbine")
            self._number_of_units = 1
        else:
            self._number_of_units = number_of_units

        # # Initialize the minimum load of the gas turbine
        # if min_load < 0.0 or min_load> 1.00:
        #     raise Exception("The minimum GT load must be between zero and one")
        # elif min_load > max_load:
        #     raise Exception("The minimum GT load must be lower than the maximum load")
        # else:
        #     self._min_load = min_load
        #
        # # Initialize the maximum load of the gas turbine
        # if max_load < 0.0 or max_load> 1.00:
        #     raise Exception("The maximum GT load must be between zero and one")
        # elif max_load < min_load:
        #     raise Exception("The maximum GT load must be higher than the minimum load")
        # else:
        #     self._max_load = max_load


    def __str__(self):
        class_info = "\nGas turbine specifications\n"
        class_info += '{:<30}{}\n'.format("GT model:",                 self.model)
        class_info += '{:<30}{}\n'.format("GT number of units:",       self.number_of_units)
        class_info += '{:<30}{:.3f}\n'.format("GT unit power (MW):",   self.unit_power/1e6)
        class_info += '{:<30}{:.3f}\n'.format("GT rated power (MW):",  self.rated_power/1e6)
        class_info += '{:<30}{:.2f}\n'.format("GT design efficiency:", self.design_efficiency*100)
        return class_info


    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model_name):
        self._model = model_name

        # Specify the power curve of each turbine model
        if model_name == 'LM2500+G4':

            # Obtained from the LM2500+G4 map at T_in=15C and p_in=1atm (Excel file)
            self._unit_power = 32.245e6
            self._load_data = [0.0000, 0.0880, 0.1035, 0.1190, 0.1344, 0.1499, 0.1653, 0.1808, 0.1962, 0.2117, 0.2271,
                               0.2425, 0.2580, 0.2734, 0.2888, 0.3042, 0.3196, 0.3351, 0.3505, 0.3659, 0.3813, 0.3967,
                               0.4121, 0.4275, 0.4429, 0.4582, 0.4736, 0.4890, 0.5044, 0.5197, 0.5351, 0.5504, 0.5658,
                               0.5811, 0.5965, 0.6118, 0.6272, 0.6425, 0.6579, 0.6732, 0.6885, 0.7038, 0.7192, 0.7345,
                               0.7498, 0.7651, 0.7804, 0.7957, 0.8110, 0.8263, 0.8416, 0.8568, 0.8721, 0.8874, 0.9027,
                               0.9179, 0.9332, 0.9485, 0.9637, 0.9790, 0.9942, 0.9999, 1.000]
            self._eff_data = [0.1297, 0.1297, 0.1424, 0.1543, 0.1679, 0.1809, 0.1934, 0.2055, 0.2171, 0.2283, 0.2385,
                              0.2456, 0.2522, 0.2315, 0.2370, 0.2422, 0.2471, 0.2519, 0.2564, 0.2607, 0.2648, 0.2687,
                              0.2736, 0.2795, 0.2853, 0.2910, 0.2965, 0.3021, 0.3075, 0.3130, 0.2998, 0.3022, 0.3046,
                              0.3068, 0.3089, 0.3110, 0.3151, 0.3194, 0.3236, 0.3277, 0.3317, 0.3356, 0.3396, 0.3435,
                              0.3474, 0.3512, 0.3550, 0.3587, 0.3624, 0.3660, 0.3696, 0.3732, 0.3759, 0.3770, 0.3782,
                              0.3791, 0.3801, 0.3810, 0.3818, 0.3823, 0.3828, 0.3830, 0.3830]

            # Obtained from the LM2500+G4 and WHRU map
            self._hload_data = [0.000, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700,
                                0.750, 0.800, 0.850, 0.900, 0.950, 1.000]
            self._heat_data  = [0.000, 11.00, 12.00, 14.00, 15.00, 16.00, 17.00, 17.50, 18.00, 19.00, 19.50, 20.00,
                                20.20, 20.40, 20.60, 20.80, 21.00, 22.00]

        elif model_name == 'LM6000-PF':

            #  Obtained from the LM6000-PF map at T_in=15C and p_in=1atm (Excel file)
            self._unit_power = 41.945e6
            self._load_data = [0.0000, 0.0902, 0.1020, 0.1139, 0.1257, 0.1376, 0.1495, 0.1613, 0.1732, 0.1850, 0.1969,
                               0.2087, 0.2206, 0.2324, 0.2443, 0.2561, 0.2680, 0.2798, 0.2917, 0.3035, 0.3154, 0.3272,
                               0.3390, 0.3509, 0.3627, 0.3746, 0.3864, 0.3983, 0.4101, 0.4219, 0.4338, 0.4456, 0.4574,
                               0.4693, 0.4811, 0.4929, 0.5048, 0.5166, 0.5284, 0.5403, 0.5521, 0.5639, 0.5757, 0.5876,
                               0.5994, 0.6112, 0.6230, 0.6348, 0.6466, 0.6585, 0.6703, 0.6821, 0.6939, 0.7057, 0.7175,
                               0.7293, 0.7411, 0.7529, 0.7647, 0.7765, 0.7883, 0.8001, 0.8119, 0.8237, 0.8355, 0.8473,
                               0.8591, 0.8709, 0.8827, 0.8945, 0.9063, 0.9181, 0.9298, 0.9416, 0.9534, 0.9652, 0.9769,
                               0.9887, 0.9999, 1.0000]
            self._eff_data = [0.1095, 0.1095, 0.1200, 0.1299, 0.1393, 0.1489, 0.1589, 0.1687, 0.1782, 0.1876, 0.1968,
                              0.1924, 0.1985, 0.2045, 0.2103, 0.2159, 0.2213, 0.2265, 0.2318, 0.2366, 0.2411, 0.2454,
                              0.2512, 0.2555, 0.2597, 0.2636, 0.2672, 0.2708, 0.2741, 0.2775, 0.2824, 0.2875, 0.2924,
                              0.2972, 0.3020, 0.3068, 0.3116, 0.3163, 0.3210, 0.3256, 0.3300, 0.3345, 0.3390, 0.3434,
                              0.3470, 0.3500, 0.3530, 0.3559, 0.3587, 0.3615, 0.3642, 0.3669, 0.3696, 0.3722, 0.3748,
                              0.3768, 0.3785, 0.3801, 0.3816, 0.3832, 0.3847, 0.3863, 0.3879, 0.3895, 0.3910, 0.3926,
                              0.3942, 0.3964, 0.3992, 0.4019, 0.4036, 0.4046, 0.4057, 0.4067, 0.4073, 0.4089, 0.4099,
                              0.4109, 0.4116, 0.4116]

            # Obtained from the LM6000-PF and WHRU map
            self._hload_data = [0.000, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700,
                                0.750, 0.800, 0.850, 0.900, 0.950, 1.000]
            self._heat_data  = [0.000, 12.00, 15.00, 16.00, 16.50, 17.00, 17.50, 18.00, 18.20, 18.40, 18.60, 18.80,
                                19.00, 19.20, 19.40, 19.60, 19.80, 20.00]

        else:
            raise Exception("Invalid gas turbine model\nValid options: 'LM2500+G4', 'LM6000-PF'")

        # Create a function to interpolate the efficiency
        self._load_data, self._eff_data = np.asarray(self._load_data), np.asarray(self._eff_data)
        self._efficiency_interpolant = interp1d(self._load_data, self._eff_data, kind='linear')

        # Create a function to interpolate the heat load
        self._load_hdata, self._heat_data = np.asarray(self._hload_data), np.asarray(self._heat_data)*1e6


    @property
    def number_of_units(self):
        return self._number_of_units

    @number_of_units.setter
    def number_of_units(self, number_of_units):
        self._number_of_units = number_of_units

    @property
    def unit_power(self):
        return self._unit_power

    @property
    def rated_power(self):
        return self._unit_power * self.number_of_units

    # @property
    # def min_power(self):
    #     return self._min_load * self.rated_power
    #
    # @property
    # def max_power(self):
    #     return self._max_load * self.rated_power

    @property
    def design_efficiency(self):
        return self._eff_data[-1]


    # def compute_power_output(self, power_demand):
    #     # Compute power output as a piece-wise function of the power demand
    #     power = np.where(np.logical_and(power_demand > 0, power_demand < self.min_power), self.min_power, 0) + \
    #             np.where(np.logical_and(power_demand >= self.min_power, power_demand < self.max_power), power_demand,0) + \
    #             np.where(np.logical_and(power_demand >= self.max_power, True), self.max_power, 0)
    #     return power  # Units: W

    # def compute_heat_output(self, power_output):
    #     # Compute heat output for a given power output
    #     efficiency = self.compute_efficiency(power_output)
    #     heat_output = power_output * (1-efficiency) / efficiency
    #     return heat_output  # Units: W

    def compute_minimum_power_output(self, heat_demand):
        # Compute minimum power output required to satisfy heat demand
        if self.compute_heat_from_power(self.rated_power) < heat_demand:
            raise Exception("The current GT configuration cannot satisfy the heat demand")
        else:
            heat_error = lambda power: heat_demand - self.compute_heat_from_power(power)
            solution = root_scalar(heat_error, method="brentq", bracket=[0.00, self.rated_power])
            minimum_power = solution.root
        return minimum_power  # Units: W

    def compute_heat_from_power(self, power_output):
        heat_interpolant = interp1d(self._hload_data, self._heat_data, kind='linear')
        heat_output = heat_interpolant(power_output/self.rated_power)
        return heat_output  # Units: W

    def compute_power_from_heat(self, heat_output):
        load_interpolant = interp1d(self._heat_data, self._hload_data, kind='linear')
        power_output = load_interpolant(heat_output)*self.rated_power
        return power_output

    def compute_efficiency(self, power_output):
        # Compute the gas turbine efficiency from the performance map
        if np.any(power_output < 0) or np.any(power_output > self.rated_power):
            raise Exception("The power output is outside the range of operation")
        load_fraction = power_output / self.rated_power
        efficiency = self._efficiency_interpolant(load_fraction)
        return efficiency  # Units: fraction of unity

    def compute_fuel_mass_flow_rate(self, power_output, fuel):
        # Compute the mass flow rate of fuel for a given power demand and fuel composition
        efficiency = self.compute_efficiency(power_output)
        mass_flow_fuel = compute_fuel_consumption(power_output, efficiency, fuel.x, fuel.MW, fuel.LHV)
        return mass_flow_fuel  # Units: kg/s

    def compute_hydrogen_mass_flow_rate(self, power_output, fuel):
        # Compute the mass flow rate of fuel for a given power demand and fuel composition
        efficiency = self.compute_efficiency(power_output)
        mass_flow_fuel = compute_fuel_consumption(power_output, efficiency, fuel.x, fuel.MW, fuel.LHV)
        mass_flow_hydrogen = mass_flow_fuel * fuel.y[np.where(fuel.components=="H2")]
        if mass_flow_hydrogen.size == 0:
            mass_flow_hydrogen = 0.0
        return mass_flow_hydrogen  # Units: kg/s

    def compute_carbon_dioxide_emissions(self, power_output, fuel):
        # Compute the mass flow rate of carbon dioxide for a given power demand and fuel composition
        efficiency = self.compute_efficiency(power_output)
        mass_flow_CO2 = compute_carbon_dioxide_emissions(power_output, efficiency, fuel.x, fuel.MW, fuel.LHV, fuel.CR)
        return mass_flow_CO2  # Units: kg/s

    def compute_specific_carbon_dioxide_emissions(self, power_output, fuel):
        # Compute the mass flow rate of carbon dioxide for a given power demand and fuel composition
        efficiency = self.compute_efficiency(power_output)
        specific_CO2 = compute_specific_carbon_dioxide_emissions(efficiency, fuel.x, fuel.MW, fuel.LHV, fuel.CR)
        return specific_CO2  # Units: g/kWh


    # def plot_power_curve(self, fig=None, ax=None):
    #     # Plot the power output as a function of the power demand
    #     if fig is None:
    #         fig = plt.figure(figsize=(5.4, 5))
    #         ax = fig.add_subplot(111)
    #     fontsize = 14
    #     ax.set_xlabel('Power demand (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
    #     ax.set_ylabel('Power output (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
    #     # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
    #     # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
    #     for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
    #     for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
    #     power_demand = np.linspace(0, self.rated_power, 200)
    #     power_output = self.compute_power_output(power_demand)
    #     ax.plot(power_demand / 1e6, power_output / 1e6, linewidth=1.25, linestyle='-', color='k',
    #             marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
    #     fig.tight_layout()
    #     return fig, ax

    def plot_efficiency_curve(self, fig=None, ax=None):
        # Plot the gas turbine efficiency against the fraction of load
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Gas turbine load (%)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Gas turbine efficiency (%)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        load = np.linspace(0, 1, 200)
        efficiency = self._efficiency_interpolant(load)
        ax.plot(load * 100, efficiency * 100, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        ax.plot(self._load_data * 100, self._eff_data * 100, linewidth=1.25, linestyle=' ', color='k',
                marker='o', markersize=3.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax

    def plot_carbon_dioxide_emissions_curve(self, fig=None, ax=None, fuel=natural_gas):
        # Plot the gas turbine carbon dioxide emissions as a function of the power demand
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Power output (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('CO$_2$ emissions (kg/s)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        power_output = np.linspace(0, self.rated_power, 200)
        CO2_emissions = self.compute_carbon_dioxide_emissions(power_output, fuel)
        ax.plot(power_output / 1e6, CO2_emissions, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w',
                label=fuel.name)
        ax.legend(ncol=1, loc='best', fontsize=fontsize - 2, edgecolor='k', framealpha=1.0)
        fig.tight_layout()
        return fig, ax

    def plot_specific_carbon_dioxide_emissions_curve(self, fig=None, ax=None, fuel=natural_gas):
        # Plot the gas turbine carbon dioxide emissions as a function of the power demand
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Power output (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Specific CO$_2$ emissions (g/kWh)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        power_output = np.linspace(0, self.rated_power, 200)
        CO2_emissions = self.compute_specific_carbon_dioxide_emissions(power_output, fuel)
        ax.plot(power_output / 1e6, CO2_emissions, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w',
                label=fuel.name)
        ax.legend(ncol=1, loc='best', fontsize=fontsize - 2, edgecolor='k', framealpha=1.0)
        fig.tight_layout()
        return fig, ax


    def plot_heat_and_power_curve(self, fig=None, ax=None, fuel=natural_gas):
        # Plot the gas turbine carbon dioxide emissions as a function of the power demand
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Gas turbine load (%)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Gas turbine heat and power (MW)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        GT_load = np.linspace(0.00, 1.00, 250)
        GT_power = GT_load*self.rated_power
        GT_heat = self.compute_heat_from_power(GT_power)
        GT_power_bis = self.compute_power_from_heat((GT_heat))
        ax.plot(GT_load*100, GT_power/1e6, linewidth=1.25, linestyle='-', color='k', marker=' ', markersize=4.5,
                markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label="Power output")
        ax.plot(GT_load*100, GT_power_bis/1e6, linewidth=1.00, linestyle='--', color='g', marker=' ', markersize=4.5,
                markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label="Power output")
        ax.plot(GT_load*100, GT_heat/1e6, linewidth=1.25, linestyle='-', color='r', marker=' ', markersize=4.5,
                markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label="Heat output")
        ax.legend(ncol=1, loc='best', fontsize=fontsize - 2, edgecolor='k', framealpha=1.0)
        fig.tight_layout()
        return fig, ax


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


class EL:
    """
    The EL object defines and stores parameters of the electrolizer system
    """

    def __init__(self, model=None, efficiency_coefficients=None, rated_power=None):

        # Initialize electrolizer polynomial coefficients
        # The coefficients can only be specified upon class initialization
        self._efficiency_coefficients = np.atleast_1d(efficiency_coefficients)

        # Initialize electrolizer model
        if model is None:
            print("The electrolizer model was not specified. Using model='NEL_HYDROGEN'")
            self.model = "NEL_HYDROGEN"
        else:
            self.model = model

        # Initialize wind farm rated power
        if rated_power is None:
            print("The rated power of the electrolizer system was not specified. Using one single electrolizer stack")
            self.rated_power = self._unit_power
        else:
            self.rated_power = rated_power


    def __str__(self):
        class_info = '\nElectrolizer system specifications\n'
        class_info += '{:<30}{}\n'.format("EL model:", self.model)
        class_info += '{:<30}{}\n'.format("EL number of units:", self.number_of_units)
        class_info += '{:<30}{:.3f}\n'.format("EL unit power (MW):", self.unit_power / 1e6)
        class_info += '{:<30}{:.3f}\n'.format("EL rated power (MW):", self.rated_power / 1e6)
        class_info += '{:<30}{:.3f}\n'.format("EL max efficiency (%):", self.maximum_efficiency*100)
        class_info += '{:<30}{:.5f}\n'.format("EL max performance (kg/MJ):", self.maximum_performance*1e6)
        return class_info


    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model_name):
        # Specify the electrolyzer model
        self._model = model_name

        # Specify the efficiency curve of the electrolizer stack
        if model_name == 'NEL_HYDROGEN':
            self._unit_power = 500e3
            self._efficiency_coefficients = 6.82299e-03, -8.29304e-09, 1.59861e-14, -1.41175e-20
            self._efficiency_coefficients = [coeff / 1e6 for coeff in self._efficiency_coefficients] # Convert units to kgH2/J
            self._efficiency_coefficients = self._efficiency_coefficients * hydrogen.LHV
            self._evaluate_stack_efficiency = \
                lambda stack_power: sum([coeff * stack_power ** i for i, coeff in enumerate(self._efficiency_coefficients)])

        elif model_name == 'POLYNOMIAL_EFFICIENCY':
            if self._efficiency_coefficients[0] is None:
                raise Exception("The coefficients must be specified when model='POLYNOMIAL_EFFICIENCY'")
            self._unit_power = 500e3
            self._efficiency_coefficients = self._efficiency_coefficients
            self._evaluate_stack_efficiency = \
                lambda stack_power: sum([coeff * stack_power ** i for i, coeff in enumerate(self._efficiency_coefficients)])

        else:
            raise Exception("Invalid electrolizer model\nValid options: 'NEL_HYDROGEN', 'POLYNOMIAL_EFFICIENCY'")

    @property
    def unit_power(self):
        return self._unit_power

    @unit_power.setter
    def unit_power(self, unit_power):
        if unit_power >= 0:
            self._unit_power = unit_power
        else:
            raise ValueError("Unit power must be positive")

    @property
    def rated_power(self):
        return self._rated_power

    @rated_power.setter
    def rated_power(self, rated_power):
        if rated_power >= 0:
            self._rated_power = rated_power
        else:
            raise ValueError("Rated power must be positive")

    @property
    def number_of_units(self):
        return self.rated_power / self.unit_power

    @property
    def maximum_efficiency(self):
        power_input = np.linspace(0, self.rated_power, 200)
        efficiency = self.compute_hydrogen_conversion_efficiency(power_input)
        return np.max(efficiency)

    @property
    def maximum_performance(self):
        power_input = np.linspace(0, self.rated_power, 200)
        performance = self.compute_specific_hydrogen_production(power_input)
        return np.max(performance)


    def compute_hydrogen_conversion_efficiency(self, power_input): # power_input units: W
        if np.any(power_input > self.rated_power):
            raise Exception("The power input cannot be larger than the rated power")
        power_per_stack = power_input / self.number_of_units if self.rated_power != 0 else 0.0
        efficiency = self._evaluate_stack_efficiency(power_per_stack)
        return efficiency   # units: fraction of unity

    def compute_specific_hydrogen_production(self, power_input): # power_input units: W
        return self.compute_hydrogen_conversion_efficiency(power_input)/hydrogen.LHV  # units: kg/J

    def compute_hydrogen_production(self, power_input): # power_input units: W
        hydrogen_mass_flow = power_input/hydrogen.LHV*self.compute_hydrogen_conversion_efficiency(power_input)
        return hydrogen_mass_flow   # units: kg/s


    def plot_hydrogen_conversion_efficiency_curve(self, fig=None, ax=None):
        # Plot the conversion efficiency curve of the electrolizer system
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Power input (kW)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Hydrogen conversion efficiency (%)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        power_input = np.linspace(0, self.rated_power, 200)
        conversion_efficiency = self.compute_hydrogen_conversion_efficiency(power_input)
        ax.plot(power_input/1e3, conversion_efficiency*100, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax

    def plot_specific_hydrogen_production_curve(self, fig=None, ax=None):
        # Plot the performance curve of the electrolizer system
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Power input (kW)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Specific hydrogen production (kg/MJ)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        power_input = np.linspace(0, self.rated_power, 200)
        hydrogen_production = self.compute_specific_hydrogen_production(power_input)
        ax.plot(power_input/1e3, hydrogen_production*1e6, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax

    def plot_hydrogen_production_curve(self, fig=None, ax=None):
        # Plot the hydrogen production curve of the electrolizer system
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Power input (kW)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Hydrogen production (kg/s)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        power_input = np.linspace(0, self.rated_power, 200)
        hydrogen_production = self.compute_hydrogen_production(power_input)
        ax.plot(power_input/1e3, hydrogen_production, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


class FC:
    """
    The FC object defines and stores parameters of the fuel cell system
    """

    def __init__(self, model=None, efficiency_coefficients=None, rated_power=None):

        # Initialize fuel cell polynomial coefficients
        # The coefficients can only be specified upon class initialization
        self._efficiency_coefficients = np.atleast_1d(efficiency_coefficients)

        # Initialize fuel cell model
        if model is None:
            print("The fuel cell model was not specified. Using model='POWERCELL_S3'")
            self.model = "POWERCELL_S3"
        else:
            self.model = model

        # Initialize wind farm rated power
        if rated_power is None:
            print("The rated power of the fuel cell system was not specified. Using one single fuel cell stack")
            self.rated_power = self._unit_power
        else:
            self.rated_power = rated_power


    def __str__(self):
        class_info = '\nFuel cell system specifications\n'
        class_info += '{:<30}{}\n'.format("FC model:", self.model)
        class_info += '{:<30}{}\n'.format("FC number of units:", self.number_of_units)
        class_info += '{:<30}{:.3f}\n'.format("FC unit power (MW):", self.unit_power / 1e6)
        class_info += '{:<30}{:.3f}\n'.format("FC rated power (MW):", self.rated_power / 1e6)
        class_info += '{:<30}{:.3f}\n'.format("FC max efficiency (%):", self.maximum_efficiency*100)
        class_info += '{:<30}{:.2f}\n'.format("FC max performance (MJ/kg):", self.maximum_performance/1e6)
        return class_info


    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model_name):

        # Specify the performance curve of the fuel cell
        # Use model='POLYNOMIAL_EFFICIENCY' to specify the conversion efficiency curve
        self._model = model_name

        if model_name == 'POWERCELL_S3':
            # The math model for POWERCELL_S3 is a special case defined within compute_specific_energy_production()
            self._unit_power = 125e3
            A = (313800 / (2 * 96485)) * (2 / 1000)
            C0, C1, C2 = 0.003237, 4.47e-06, 7.767e-12
            self._evaluate_stack_efficiency =\
                lambda stack_power: 1e-9 + stack_power / (A * (C0 + C1 * stack_power + C2 * stack_power ** 2)) / hydrogen.LHV

        elif model_name == 'POLYNOMIAL_EFFICIENCY':
            if self._efficiency_coefficients[0] is None:
                raise Exception("The coefficients must be specified when model='POLYNOMIAL_EFFICIENCY'")
            self._unit_power = 125e3
            self._efficiency_coefficients = self._efficiency_coefficients
            self._evaluate_stack_efficiency = \
                lambda stack_power: sum([coeff * stack_power ** i for i, coeff in enumerate(self._efficiency_coefficients)])

        else:
            raise Exception("Invalid electrolizer model\nValid options: 'POWERCELL_S3', 'POLYNOMIAL_EFFICIENCY'")

    @property
    def unit_power(self):
        return self._unit_power

    @unit_power.setter
    def unit_power(self, unit_power):
        if unit_power >= 0:
            self._unit_power = unit_power
        else:
            raise ValueError("Unit power must be positive")

    @property
    def rated_power(self):
        return self._rated_power

    @rated_power.setter
    def rated_power(self, rated_power):
        if rated_power >= 0:
            self._rated_power = rated_power
        else:
            raise ValueError("Rated power must be positive")

    @property
    def number_of_units(self):
        return self.rated_power / self.unit_power

    @property
    def maximum_efficiency(self):
        power_input = np.linspace(0, self.rated_power, 200)
        efficiency = self.compute_hydrogen_conversion_efficiency(power_input)
        return np.max(efficiency)

    @property
    def maximum_performance(self):
        power_input = np.linspace(0, self.rated_power, 200)
        performance = self.compute_specific_energy_production(power_input)
        return np.max(performance)


    def compute_hydrogen_conversion_efficiency(self, power_output): # power_output units: W
        if np.any(power_output > self.rated_power):
            raise Exception("The power output cannot be larger than the rated power")
        power_per_stack = power_output / self.number_of_units if self.rated_power != 0 else 0.0
        efficiency = self._evaluate_stack_efficiency(power_per_stack)
        return efficiency   # units: fraction of unity

    def compute_specific_energy_production(self, power_output): # power_output units: W
        return self.compute_hydrogen_conversion_efficiency(power_output)*hydrogen.LHV   # units: J/kg

    def compute_hydrogen_consumption(self, power_output): # power_output units: W
        hydrogen_mass_flow = power_output / hydrogen.LHV / self.compute_hydrogen_conversion_efficiency(power_output)
        return hydrogen_mass_flow   # units: kg/s


    def plot_hydrogen_conversion_efficiency_curve(self, fig=None, ax=None):
        # Plot the conversion efficiency curve of the fuel cell system
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Power output (kW)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Hydrogen conversion efficiency (%)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        power_output = np.linspace(0, self.rated_power, 200)
        conversion_efficiency = self.compute_hydrogen_conversion_efficiency(power_output)
        ax.plot(power_output/1e3, conversion_efficiency*100, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax

    def plot_hydrogen_consumption_curve(self, fig=None, ax=None):
        # Plot the hydrogen consumption curve of the fuel cell system
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Power output (kW)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Hydrogen consumption (kg/s)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        power_output = np.linspace(0, self.rated_power, 200)
        hydrogen_consumption = self.compute_hydrogen_consumption(power_output)
        ax.plot(power_output/1e3, hydrogen_consumption, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax

    def plot_specific_electricity_production_curve(self, fig=None, ax=None):
        # Plot the performance curve of the fuel cell system
        if fig is None:
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111)
        fontsize = 14
        ax.set_xlabel('Power output (kW)', fontsize=fontsize, color='k', labelpad=fontsize)
        ax.set_ylabel('Specific electricity production (MJ/kg)', fontsize=fontsize, color='k', labelpad=fontsize)
        # ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
        for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
        power_output = np.linspace(0, self.rated_power, 200)
        hydrogen_consumption = self.compute_specific_energy_production(power_output)
        ax.plot(power_output/1e3, hydrogen_consumption/1e6, linewidth=1.25, linestyle='-', color='k',
                marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=None)
        fig.tight_layout()
        return fig, ax



## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##

