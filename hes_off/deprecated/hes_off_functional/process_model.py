# Import packages
import numpy as np
import numba as nb
from scipy.io import loadmat
from importlib_resources import files

# Compile time constants
hydrogen_LHV = 119.96e6
MW_CO2 = 44.01/1e3


## ------------------------------------------------------------------------------------------------------------------ ##
## Wind turbine functions
## ------------------------------------------------------------------------------------------------------------------ ##

# HYWIND turbine data based on based on https://www.uib.no/sites/w3.uib.no/files/attachments/hywind_energy_lab.pdf (page 25)
HYWIND_unit_power = 6e6
HYWIND_data_power = np.asarray(((
                    0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00,
                    11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 25.00, 26.00, 50.00),
                    (0.0000, 0.0000, 0.0000, 0.0000, 0.0290, 0.0725, 0.1304, 0.2101, 0.3261, 0.4638,
                    0.6232, 0.7754, 0.8913, 0.9565, 0.9855, 1.0000, 1.0000, 1.0000, 0.0000, 0.0000)))

# NREL turbine data based on report NREL/TP-500-38060
NREL_unit_power = 5e6
NREL_data_power = np.asarray(((
                    0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00,
                    11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 25.00, 30.00, 100.00),
                    (0, 0, 0, 170, 391, 731, 1173, 1752, 2534, 3452, 4558, 5000, 5000, 5000,
                    5000, 5000, 5000, 5000, 0, 0)))


# @nb.jit(nopython=True, cache=True)
def compute_WT_power_output(model, rated_power, hub_height, ref_height, wind_speed):

    if model == 'HYWIND':
        unit_power = HYWIND_unit_power
        speed_data = HYWIND_data_power[0,:]
        power_data = HYWIND_data_power[1,:]

    elif model == 'NREL':
        unit_power = NREL_unit_power
        speed_data = NREL_data_power[0,:]
        power_data = NREL_data_power[1,:]
    else:
        raise Exception("Invalid wind turbine model\nValid options: 'HYWIND', 'NREL'")

    # Normalize the power data
    normalized_power_data = power_data / np.max(power_data)
    power_data = normalized_power_data * rated_power

    # Compute the wind speed factor
    wind_speed_factor = (hub_height / ref_height) ** (1 / 7)

    # Compute the wind turbine power output
    power_output = np.interp(wind_speed*wind_speed_factor, speed_data, power_data)

    return power_output


# @nb.jit(nopython=True, cache=True)
def compute_WT_energy_output(model, rated_power, hub_height, ref_height, wind_speed, time_series):
    power_timeseries =  compute_WT_power_output(model, rated_power, hub_height, ref_height, wind_speed)
    integrated_energy = np.trapz(power_timeseries, time_series * 3600)
    return integrated_energy


# @nb.jit(nopython=True, cache=True)
def compute_WT_capacity_factor(model, rated_power, hub_height, ref_height, wind_speed, time_series):
    integrated_energy = compute_WT_energy_output(model, rated_power, hub_height, ref_height, wind_speed, time_series)
    return integrated_energy / (rated_power * 3600 * 24 * 365)


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


## ------------------------------------------------------------------------------------------------------------------ ##
## Electrolyzer functions
## ------------------------------------------------------------------------------------------------------------------ ##

# @nb.jit(nopython=True, cache=True)
def compute_EL_efficiency(model, efficiency_coefficients, rated_power, power_input): # power_output units: W

    # Convert arguments to numpy arrays
    power_input = np.atleast_1d(np.asarray(power_input))

    # Check power input bounds
    if np.any(power_input < 0.0) or np.any(power_input > rated_power):
        raise Exception("The power input is outside the range of operation")

    # Define the fuel cell model
    if model == 'NEL_HYDROGEN':
        # The math model for POWERCELL_S3 is a special case defined within compute_specific_energy_production()
        unit_power = 500e3
        number_of_units = rated_power / unit_power
        stack_power = power_input / number_of_units if rated_power != 0 else 0.0*power_input
        efficiency_coefficients = (6.82299e-03, -8.29304e-09, 1.59861e-14, -1.41175e-20)
        efficiency = 0.0*power_input
        for i, coeff in enumerate(efficiency_coefficients):
            efficiency += (coeff * stack_power ** i) * hydrogen_LHV / 1e6

    elif model == 'POLYNOMIAL_EFFICIENCY':
        unit_power = 500e3
        number_of_units = rated_power / unit_power
        stack_power = power_input / number_of_units if rated_power != 0 else 0.0*power_input
        efficiency = 0.0*power_input
        for i, coeff in enumerate(efficiency_coefficients):
            efficiency += (coeff * stack_power ** i) * hydrogen_LHV / 1e6
    else:
        raise Exception("Invalid electrolizer model\nValid options: 'NEL_HYDROGEN', 'POLYNOMIAL_EFFICIENCY'")

    return efficiency  # units: fraction of unity


# @nb.jit(nopython=True, cache=True)
def compute_EL_specific_performance(model, efficiency_coefficients, rated_power, power_input): # power_output units: W
    efficiency = compute_EL_efficiency(model, efficiency_coefficients, rated_power, power_input)
    return efficiency / hydrogen_LHV


# @nb.jit(nopython=True, cache=True)
def compute_EL_hydrogen_production(model, efficiency_coefficients, rated_power, power_input): # power_output units: W
    efficiency = compute_EL_efficiency(model, efficiency_coefficients, rated_power, power_input)
    return power_input / hydrogen_LHV * efficiency



## ------------------------------------------------------------------------------------------------------------------ ##
## Fuel cell functions
## ------------------------------------------------------------------------------------------------------------------ ##

# @nb.jit(nopython=True, cache=True)
def compute_FC_efficiency(model, efficiency_coefficients, rated_power, power_output): # power_output units: W

    # Convert arguments to numpy arrays
    power_output = np.atleast_1d(np.asarray(power_output))

    # Check power output bounds
    if np.any(power_output < 0.0) or np.any(power_output > rated_power):
        raise Exception("The power input is outside the range of operation")

    # Define the fuel cell model
    if model == 'POWERCELL_S3':
        # The math model for POWERCELL_S3 is a special case defined within compute_specific_energy_production()
        unit_power = 125e3
        number_of_units = rated_power / unit_power
        stack_power = power_output / number_of_units if rated_power != 0 else 0.0*power_output
        A = (313800 / (2 * 96485)) * (2 / 1000)
        C0, C1, C2 = (0.003237, 4.47e-06, 7.767e-12)
        efficiency = 1e-9 + stack_power / (A * (C0 + C1 * stack_power + C2 * stack_power ** 2)) / hydrogen_LHV

    elif model == 'POLYNOMIAL_EFFICIENCY':
        unit_power = 125e3
        number_of_units = rated_power / unit_power
        stack_power = power_output / number_of_units if rated_power != 0 else 0.0*power_output
        efficiency = 0.0*power_output
        for i, coeff in enumerate(efficiency_coefficients):
            efficiency += (coeff * stack_power ** i) * hydrogen_LHV / 1e6
    else:
        raise Exception("Invalid electrolizer model\nValid options: 'POWERCELL_S3', 'POLYNOMIAL_EFFICIENCY'")

    return efficiency  # units: fraction of unity


# @nb.jit(nopython=True, cache=True)
def compute_FC_specific_performance(model, efficiency_coefficients, rated_power, power_output): # power_output units: W
    efficiency = compute_FC_efficiency(model, efficiency_coefficients, rated_power, power_output)
    return  efficiency * hydrogen_LHV   # units: J/kg


# @nb.jit(nopython=True, cache=True)
def compute_FC_hydrogen_consumption(model, efficiency_coefficients, rated_power, power_output): # power_output units: W
    efficiency = compute_FC_efficiency(model, efficiency_coefficients, rated_power, power_output)
    return power_output / hydrogen_LHV / efficiency   # units: kg/s



## ------------------------------------------------------------------------------------------------------------------ ##
## Gas turbine functions
## ------------------------------------------------------------------------------------------------------------------ ##

# Obtained from the LM2500+G4 map at T_in=15C and p_in=1atm (Excel file)
LM2500_unit_power = 32.245e6
LM2500_unit_heat = 22.000e6
LM2500_data_efficiency = np.asarray(((
                            0.0000, 0.0880, 0.1035, 0.1190, 0.1344, 0.1499, 0.1653, 0.1808, 0.1962, 0.2117, 0.2271,
                            0.2425, 0.2580, 0.2734, 0.2888, 0.3042, 0.3196, 0.3351, 0.3505, 0.3659, 0.3813, 0.3967,
                            0.4121, 0.4275, 0.4429, 0.4582, 0.4736, 0.4890, 0.5044, 0.5197, 0.5351, 0.5504, 0.5658,
                            0.5811, 0.5965, 0.6118, 0.6272, 0.6425, 0.6579, 0.6732, 0.6885, 0.7038, 0.7192, 0.7345,
                            0.7498, 0.7651, 0.7804, 0.7957, 0.8110, 0.8263, 0.8416, 0.8568, 0.8721, 0.8874, 0.9027,
                            0.9179, 0.9332, 0.9485, 0.9637, 0.9790, 0.9942, 0.9999, 1.000),
                            (0.1297, 0.1297, 0.1424, 0.1543, 0.1679, 0.1809, 0.1934, 0.2055, 0.2171, 0.2283, 0.2385,
                            0.2456, 0.2522, 0.2315, 0.2370, 0.2422, 0.2471, 0.2519, 0.2564, 0.2607, 0.2648, 0.2687,
                            0.2736, 0.2795, 0.2853, 0.2910, 0.2965, 0.3021, 0.3075, 0.3130, 0.2998, 0.3022, 0.3046,
                            0.3068, 0.3089, 0.3110, 0.3151, 0.3194, 0.3236, 0.3277, 0.3317, 0.3356, 0.3396, 0.3435,
                            0.3474, 0.3512, 0.3550, 0.3587, 0.3624, 0.3660, 0.3696, 0.3732, 0.3759, 0.3770, 0.3782,
                            0.3791, 0.3801, 0.3810, 0.3818, 0.3823, 0.3828, 0.3830, 0.3830)))

# Obtained from the LM2500+G4 + WHRU map
LM2500_data_heat = np.asarray(((
                            0.000, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550,
                            0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000),
                            (0.000, 11.00, 12.00, 14.00, 15.00, 16.00, 17.00, 17.50, 18.00,
                            19.00, 19.50, 20.00, 20.20, 20.40, 20.60, 20.80, 21.00, 22.00)))
LM2500_data_heat[1,:] *= 1e6


# Obtained from the LM6000-PF map at T_in=15C and p_in=1atm (Excel file)
LM6000_unit_power = 41.945e6
LM6000_unit_heat = 20.000e6
LM6000_data_efficiency = np.asarray(((
                            0.0000, 0.0902, 0.1020, 0.1139, 0.1257, 0.1376, 0.1495, 0.1613, 0.1732, 0.1850, 0.1969,
                            0.2087, 0.2206, 0.2324, 0.2443, 0.2561, 0.2680, 0.2798, 0.2917, 0.3035, 0.3154, 0.3272,
                            0.3390, 0.3509, 0.3627, 0.3746, 0.3864, 0.3983, 0.4101, 0.4219, 0.4338, 0.4456, 0.4574,
                            0.4693, 0.4811, 0.4929, 0.5048, 0.5166, 0.5284, 0.5403, 0.5521, 0.5639, 0.5757, 0.5876,
                            0.5994, 0.6112, 0.6230, 0.6348, 0.6466, 0.6585, 0.6703, 0.6821, 0.6939, 0.7057, 0.7175,
                            0.7293, 0.7411, 0.7529, 0.7647, 0.7765, 0.7883, 0.8001, 0.8119, 0.8237, 0.8355, 0.8473,
                            0.8591, 0.8709, 0.8827, 0.8945, 0.9063, 0.9181, 0.9298, 0.9416, 0.9534, 0.9652, 0.9769,
                            0.9887, 0.9999, 1.0000),
                            (0.1095, 0.1095, 0.1200, 0.1299, 0.1393, 0.1489, 0.1589, 0.1687, 0.1782, 0.1876, 0.1968,
                            0.1924, 0.1985, 0.2045, 0.2103, 0.2159, 0.2213, 0.2265, 0.2318, 0.2366, 0.2411, 0.2454,
                            0.2512, 0.2555, 0.2597, 0.2636, 0.2672, 0.2708, 0.2741, 0.2775, 0.2824, 0.2875, 0.2924,
                            0.2972, 0.3020, 0.3068, 0.3116, 0.3163, 0.3210, 0.3256, 0.3300, 0.3345, 0.3390, 0.3434,
                            0.3470, 0.3500, 0.3530, 0.3559, 0.3587, 0.3615, 0.3642, 0.3669, 0.3696, 0.3722, 0.3748,
                            0.3768, 0.3785, 0.3801, 0.3816, 0.3832, 0.3847, 0.3863, 0.3879, 0.3895, 0.3910, 0.3926,
                            0.3942, 0.3964, 0.3992, 0.4019, 0.4036, 0.4046, 0.4057, 0.4067, 0.4073, 0.4089, 0.4099,
                            0.4109, 0.4116, 0.4116)))

# Obtained from the LM6000-PF + WHRU map
LM6000_data_heat = np.asarray(((
                            0.000, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550,
                            0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000),
                            (0.000, 12.00, 15.00, 16.00, 16.00, 17.00, 17.50, 18.00, 18.20,
                            18.40, 18.60, 18.80, 19.00, 19.20, 19.40, 19.60, 19.80, 20.00)))
LM6000_data_heat[1,:] *= 1e6


# @nb.jit(nopython=True, cache=True)
def compute_GT_efficiency(model, number_of_units, power_output):

    # Convert arguments to numpy arrays
    power_output = np.atleast_1d(np.asarray(power_output))

    # Specify the power curve of each turbine model
    if model == 'LM2500+G4':
        unit_power = LM2500_unit_power
        rated_power = unit_power*number_of_units
        load_data = LM2500_data_efficiency[0,:]
        eff_data = LM2500_data_efficiency[1,:]

    elif model == 'LM6000-PF':
        unit_power = LM6000_unit_power
        rated_power = unit_power*number_of_units
        load_data = LM6000_data_efficiency[0,:]
        eff_data = LM6000_data_efficiency[1,:]

    else:
        raise Exception("Invalid gas turbine model\nValid options: 'LM2500+G4', 'LM6000-PF'")

    # Check power output bounds
    if np.any(power_output < 0.0) or np.any(power_output > rated_power):
        raise Exception("The power output is outside the range of operation")

    # Compute the gas turbine efficiency
    load_fraction = power_output / rated_power
    efficiency = np.interp(load_fraction, load_data, eff_data)

    return efficiency


# @nb.jit(nopython=True, cache=True)
def compute_GT_power_from_heat(model, number_of_units, heat_output):

    # Convert arguments to numpy arrays
    heat_output = np.atleast_1d(np.asarray(heat_output))

    # Specify the heat curve of each turbine model
    if model == 'LM2500+G4':
        unit_power = LM2500_unit_power
        rated_power = unit_power*number_of_units
        unit_heat = LM2500_unit_heat
        rated_heat = unit_heat*number_of_units
        load_data = LM2500_data_heat[0,:]
        heat_data = LM2500_data_heat[1,:]*number_of_units

    elif model == 'LM6000-PF':
        unit_power = LM6000_unit_power
        rated_power = unit_power*number_of_units
        unit_heat = LM6000_unit_heat
        rated_heat = unit_heat*number_of_units
        load_data = LM6000_data_heat[0,:]
        heat_data = LM6000_data_heat[1,:]*number_of_units

    else:
        raise Exception("Invalid gas turbine model\nValid options: 'LM2500+G4', 'LM6000-PF'")

    # Check power output bounds
    if np.any(heat_output < 0.0) or np.any(heat_output > rated_heat):
        raise Exception("The heat output is outside the range of operation")

    # Compute the gas turbine power output
    power_output = np.interp(heat_output, heat_data, load_data)*rated_power

    return power_output


# @nb.jit(nopython=True, cache=True)
def compute_GT_maximum_power(model, number_of_units):
    if model == 'LM2500+G4':
        return LM2500_unit_power*number_of_units
    elif model == 'LM6000-PF':
        return LM6000_unit_power*number_of_units
    else:
        raise Exception("Invalid gas turbine model\nValid options: 'LM2500+G4', 'LM6000-PF'")


# @nb.jit(nopython=True, cache=True)
def compute_GT_fuel_consumption(power_output, conversion_efficiency, fuel):
    x, MW, LHV = fuel["x"], fuel["MW"], fuel["LHV"]
    mass_flow_fuel = np.sum(x * MW, axis=0) / np.sum(x * MW * LHV, axis=0) * (power_output / conversion_efficiency)
    return mass_flow_fuel   # Units: kg/s


# @nb.jit(nopython=True, cache=True)
def compute_GT_carbon_dioxide_emissions(power_output, conversion_efficiency, fuel):
    x, MW, LHV, CR = fuel["x"], fuel["MW"], fuel["LHV"], fuel["CR"]
    mass_flow_CO2 = np.sum(CR * x, axis=0) / np.sum(x * MW * LHV, axis=0) * (power_output / conversion_efficiency) * MW_CO2
    return mass_flow_CO2   # Units: kg/s


# @nb.jit(nopython=True, cache=True)
def compute_GT_specific_carbon_dioxide_emissions(conversion_efficiency, fuel):
    mass_flow_CO2 = compute_GT_carbon_dioxide_emissions(1, conversion_efficiency, fuel)
    specific_emissions = mass_flow_CO2
    return specific_emissions * 1000 * (1000 * 3600)   # Units: g/kWh



## ------------------------------------------------------------------------------------------------------------------ ##
## Fuel mixture functions
## ------------------------------------------------------------------------------------------------------------------ ##

# @nb.jit(nopython=True, cache=True)
def convert_molar_to_mass_fraction(x, MW):
    x, MW = np.asarray(x), np.asarray(MW)
    if np.abs(1.0 - np.sum(x)) > 1e-12:
        raise Exception("The sum of the input molar fractions is not unity")
    y = (x*MW)/np.sum(x*MW)
    return y


# @nb.jit(nopython=True, cache=True)
def convert_mass_to_molar_fraction(y, MW):
    y, MW = np.asarray(y), np.asarray(MW)
    if np.abs(1.0 - np.sum(y)) > 1e-12:
        raise Exception("The sum of the input mass fractions is not unity")
    x = (y/MW)/np.sum(y/MW)
    return x


# @nb.jit(nopython=True, cache=True)
def compute_mixture_molecular_weight(x, MW):
    MW_mix = np.sum(x*MW, axis=0)
    return MW_mix   # Units: kg/mol


# @nb.jit(nopython=True, cache=True)
def create_fluid(fractions, fraction_type, MW, LHV, CR):

    # Check if the composition fractions are consistent
    if np.abs(1.0 - np.sum(fractions)) > 1e-12:
        raise Exception("The sum of the input fractions is not unity")

    # Define the fluid composition
    if fraction_type == "molar":
        x = fractions
        y = convert_molar_to_mass_fraction(x, MW)
    elif fraction_type == "mass":
        y = fractions
        x = convert_mass_to_molar_fraction(y, MW)
    else:
        raise Exception("Invalid option for 'fraction_type'. Valid options: 'molar', 'mass'")

    # Create the fluid dictionary
    fluid = {"x":          np.atleast_1d(x),
             "y":          np.atleast_1d(y),
             "CR":         np.atleast_1d(CR),
             "MW":         np.atleast_1d(MW),
             "LHV":        np.atleast_1d(LHV),
             "MW_mean":    np.atleast_1d(np.asarray(np.sum(x*MW))),
             "LHV_mean":   np.atleast_1d(np.asarray(np.sum(y*LHV)))}

    return fluid


# @nb.jit(nopython=True, cache=True)
def create_fluid_mixture(fluids, fractions, fraction_type="molar"):

    # Check if the composition fractions are consistent
    if np.abs(1.0 - np.sum(np.asarray(fractions))) > 1e-12:
        raise Exception("The sum of the input fractions is not unity")

    # Concatenate fluid properties
    x = np.empty((0,))
    y = np.empty((0,))
    MW = np.empty((0,))
    LHV = np.empty((0,))
    CR = np.empty((0,))
    for fluid, fraction in zip(fluids, fractions):
        MW = np.concatenate((MW, fluid["MW"]))
        LHV = np.concatenate((LHV, fluid["LHV"]))
        CR = np.concatenate((CR, fluid["CR"]))
        x = np.concatenate((x, fraction * fluid["x"]))
        y = np.concatenate((y, fraction * fluid["y"]))

    # Compute new molar fractions
    if fraction_type == "molar":
        y = convert_molar_to_mass_fraction(x, MW)
    elif fraction_type == "mass":
        x = convert_mass_to_molar_fraction(y, MW)
    else:
        raise Exception("Invalid option for 'fraction_type'. Valid options: 'molar', 'mass'")

    return create_fluid(fractions=x, fraction_type="molar", MW=MW, LHV=LHV, CR=CR)


def print_fluid(fluid):
    # print("{} {}\n".format("Fluid name:", fluid["name"]))
    print()
    print(("{:25}" + "{:>8}"*len(fluid["x"]) + "{:>16}").format("Components:", *range(1,len(fluid["x"])+1), "Average value"))
    print(("{:25}" + "{:>8.2f}" * len(fluid["x"]) + "{:>16.2f}").format("Molecular mass (g/mol):", *fluid["MW"]*1e3, fluid["MW_mean"][0]*1e3))
    print(("{:25}" + "{:>8.2f}" * len(fluid["x"]) + "{:>16.2f}").format("Heating value (MJ/kg):", *fluid["LHV"]/1e6, fluid["LHV_mean"][0]/1e6))
    print(("{:25}" + "{:>8.2f}" * len(fluid["x"])).format("Carbon ratio (-):", *fluid["CR"]))
    print(("{:25}" + "{:>8.4f}" * len(fluid["x"])).format("Molar composition (-):", *fluid["x"]))
    print(("{:25}" + "{:>8.4f}" * len(fluid["x"])).format("Mass composition (-):",  *fluid["y"]))



# Define the properties and composition of natural gas
natural_gas = create_fluid(
                MW=np.asarray([16.04, 30.07, 44.10, 58.12, 72.15, 86.18, 28.01, 44.01]) / 1e3,
                LHV=np.asarray([55.50, 51.90, 50.35, 49.50, 48.60, 47.7, 0.00, 0.00]) * 1e6,
                CR=np.asarray([1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 0.00, 1.00]),
                fractions=np.asarray([0.72890, 0.13601, 0.08272, 0.02750, 0.00667, 0.00091, 0.01562, 0.00167]),
                fraction_type="molar")

# Define properties and composition of methane
methane = create_fluid(
                MW=np.asarray([16.04]) / 1e3,
                LHV=np.asarray([55.50]) * 1e6,
                CR=np.asarray([1.00]),
                fractions=np.asarray([1.00]),
                fraction_type="molar")

# Define properties and composition of hydrogen
hydrogen = create_fluid(
                MW=np.asarray([2.02]) / 1e3,
                LHV=np.asarray([hydrogen_LHV]),
                CR=np.asarray([0.00]),
                fractions=np.asarray([1.00]),
                fraction_type="molar")

# # Define the properties and composition of carbon monoxide
# carbon_monoxide = create_fluid(
#                 MW=np.asarray([28.01]) / 1e3,
#                 LHV=np.asarray([10.11]) * 1e6,
#                 CR=np.asarray([1.00]),
#                 fractions=np.asarray([1.00]),
#                 fraction_type="molar")
#
# # Define the properties and composition of carbon dioxide
# carbon_dioxide = create_fluid(
#                 MW=np.asarray([MW_CO2]),
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([1.00]),
#                 fractions=np.asarray([1.00]),
#                 fraction_type="molar")
#
# # Define the properties and composition of nitrogen
# nitrogen = create_fluid(
#                 MW=np.asarray([28.01]) / 1e3,
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([0.00]),
#                 fractions=np.asarray([1.00]),
#                 fraction_type="molar")
#
# # Define the properties and composition of oxygen
# oxygen = create_fluid(
#                 MW=np.asarray([32.00]) / 1e3,
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([0.00]),
#                 fractions=np.asarray([1.00]),
#                 fraction_type="molar")
#
# # Define the properties and composition of water
# water = create_fluid(
#                 MW=np.asarray([18.02]) / 1e3,
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([0.00]),
#                 fractions=np.asarray([1.00]),
#                 fraction_type="molar")
#
# # Define the properties and composition of argon
# argon = create_fluid(
#                 MW=np.asarray([39.95]) / 1e3,
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([0.00]),
#                 fractions=np.asarray([1.00]),
#                 fraction_type="molar")
#
# # Define the properties and composition of air
# air = create_fluid(
#                 MW=np.asarray([28.01, 32.00, 44.01, 18.02, 39.95]) / 1e3,
#                 LHV=np.asarray([0.00, 0.00, 0.00, 0.00, 0.00]) * 1e6,
#                 CR=np.asarray([0.00, 0.00, 1.00, 0.00, 0.00]),
#                 fractions=np.asarray([0.7739, 0.2076, 0.0003, 0.0089, 0.0093]),
#                 fraction_type="molar")



## ------------------------------------------------------------------------------------------------------------------ ##
## End of file
## ------------------------------------------------------------------------------------------------------------------ ##
