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
import copy
import numpy as np


import numba as nb
from numba import jit
from numba.core import types
from numba.typed import Dict
from numba.experimental import jitclass

## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


float_array = types.float64[:]


@jit(nopython=True, cache=True)
def create_fluid(name, components, fractions, fraction_type, MW, LHV, CR):

    # Check if the composition fractions are consistent
    # fraction = fractions.astype(np.float64)
    # if np.abs(1.0 - np.sum(fractions)) > 1e-12:
    #     raise Exception("The sum of the input fractions is not unity")

    # Define the fluid composition
    if fraction_type == "molar":
        x = fractions
        y = convert_molar_to_mass_fraction(x, MW)
    elif fraction_type == "mass":
        y = fractions
        x = convert_mass_to_molar_fraction(y, MW)
    # else:
    #     raise Exception("Invalid option for 'fraction_type'. Valid options: 'molar', 'mass'")

    # # Create the fluid dictionary
    # fluid = {"name":       name,
    #          "components": np.atleast_1d(components),
    #          "x":          np.atleast_1d(x),
    #          "y":          np.atleast_1d(y),
    #          "CR":         np.atleast_1d(CR),
    #          "MW":         np.atleast_1d(MW),
    #          "LHV":        np.atleast_1d(LHV),
    #          "MW_mean":    np.sum(x*MW),
    #          "LHV_mean":   np.sum(y*LHV)}

    fluid = Dict.empty(
        key_type=types.unicode_type,
        value_type=float_array)
    # fluid = Dict()

    # fluid["name"] = name
    # fluid["x"] = np.atleast_1d(x).astype(np.float64)
    # fluid["y"] = np.atleast_1d(y).astype(np.float64)

    print(fluid)

    # fluid

    # return cleanup_fluid_duplicates(fluid)
    return fluid


def create_fluid_mixture(name, fluids, fractions, fraction_type="molar"):

    # Check if the composition fractions are consistent
    if np.abs(1.0 - np.sum(fractions)) > 1e-12:
        raise Exception("The sum of the input fractions is not unity")

    # Concatenate fluid properties
    MW = np.concatenate([fluid["MW"]  for fluid in fluids])
    LHV = np.concatenate([fluid["LHV"] for fluid in fluids])
    CR = np.concatenate([fluid["CR"]  for fluid in fluids])
    components = np.concatenate([fluid["components"] for fluid in fluids])

    # Compute new molar fractions
    if fraction_type == "molar":
        x = np.concatenate([fraction*fluid["x"] for fraction, fluid in zip(fractions, fluids)])
        y = convert_molar_to_mass_fraction(x, MW)
    elif fraction_type == "mass":
        y = np.concatenate([fraction*fluid["y"] for fraction, fluid in zip(fractions, fluids)])
        x = convert_mass_to_molar_fraction(y, MW)
    else:
        raise Exception("Invalid option for 'fraction_type'. Valid options: 'molar', 'mass'")

    return create_fluid(name, components, MW, LHV, CR, fraction=x, fraction_type="molar")


@jit(nopython=True, cache=True)
def cleanup_fluid_duplicates(fluid):

    # Rename variables
    # fluid = copy.deepcopy(fluid)
    components = fluid["components"]
    x = fluid["x"]
    y = fluid["y"]
    MW = fluid["MW"]
    LHV = fluid["LHV"]
    CR = fluid["CR"]

    # Iterate over the component names
    for name in components:

        # Find the entries that are duplicated
        repeated_indices = [i for i, name_bis in enumerate(components) if name == name_bis]

        # Eliminate repeated information
        components = np.delete(components, repeated_indices[1::])
        MW  = np.delete(MW, repeated_indices[1::])
        LHV = np.delete(LHV, repeated_indices[1::])
        CR  = np.delete(CR, repeated_indices[1::])

        # Adjust the molar fraction of repeated entries
        # a = x[repeated_indices[1::]]
        # print(repeated_indices)
        # print(repeated_indices[1::])
        # print(x[repeated_indices[1::]])
        # pdb.set_trace()
        x[repeated_indices[0]] = x[repeated_indices[0]] + np.sum(x[repeated_indices[1::]])
        x = np.delete(x, repeated_indices[1::])

        # Adjust the mass fraction of repeated entries
        y[repeated_indices[0]] = y[repeated_indices[0]] + np.sum(y[repeated_indices[1::]])
        y = np.delete(y, repeated_indices[1::])

    # Create the fluid dictionary
    fluid = {"name": fluid["name"],
             "components": np.atleast_1d(components),
             "x": np.atleast_1d(x),
             "y": np.atleast_1d(y),
             "CR": np.atleast_1d(CR),
             "MW": np.atleast_1d(MW),
             "LHV": np.atleast_1d(LHV),
             "MW_mean": np.sum(x*MW),
             "LHV_mean": np.sum(y*LHV)}

    return fluid


def print_fluid(fluid):
    print()
    print("hello")
    print(fluid)

    # print("{} {}\n".format("Fluid name:", fluid["name"]))
    # print(("{:30}" + "{:>8}"*len(fluid["components"]) + "{:>16}" + "\n").format("Components:", *fluid["components"], "Average value"))
    # print(("{:30}" + "{:>8.2f}" * len(fluid["components"]) + "{:>16.2f}" + "\n").format("Molecular mass (g/mol):", *fluid["MW"]*1e3, fluid["MW_mean"]*1e3))
    # print(("{:30}" + "{:>8.2f}" * len(fluid["components"]) + "{:>16.2f}" + "\n").format("Heating value (MJ/kg):", *fluid["LHV"]/1e6, fluid["LHV_mean"]/1e6))
    # print(("{:30}" + "{:>8.2f}" * len(fluid["components"]) + "\n").format("Carbon ratio (-):", *fluid["CR"]))
    # print(("{:30}" + "{:>8.4f}" * len(fluid["components"]) + "\n").format("Molar composition (-):", *fluid["x"]))
    # print(("{:30}" + "{:>8.4f}" * len(fluid["components"]) + "\n").format("Mass composition (-):",  *fluid["y"]))











@jit(nopython=True, cache=True)
def convert_molar_to_mass_fraction(x, MW):
    x, MW = np.asarray(x), np.asarray(MW)
    if np.abs(1.0 - np.sum(x)) > 1e-12:
        raise Exception("The sum of the input molar fractions is not unity")
    y = (x*MW)/np.sum(x*MW)
    return y

@jit(nopython=True, cache=True)
def convert_mass_to_molar_fraction(y, MW):
    y, MW = np.asarray(y), np.asarray(MW)
    if np.abs(1.0 - np.sum(y)) > 1e-12:
        raise Exception("The sum of the input mass fractions is not unity")
    x = (y/MW)/np.sum(y/MW)
    return x

def compute_mixture_molecular_weight(x, MW):
    MW_mix = np.sum(x*MW, axis=0)
    return MW_mix   # Units: kg/mol

def compute_fuel_consumption(power_output, conversion_efficiency, x, MW, LHV):
    mass_flow_fuel = np.sum(x * MW, axis=0) / np.sum(x * MW * LHV, axis=0) * (power_output / conversion_efficiency)
    return mass_flow_fuel   # Units: kg/s

def compute_carbon_dioxide_emissions(power_output, conversion_efficiency, x, MW, LHV, CR):
    MW_CO2 = carbon_dioxide["MW"]
    mass_flow_CO2 = np.sum(CR * x, axis=0) / np.sum(x * MW * LHV, axis=0) * (power_output / conversion_efficiency) * MW_CO2
    return mass_flow_CO2   # Units: kg/s

def compute_specific_carbon_dioxide_emissions(conversion_efficiency, x, MW, LHV, CR):
    mass_flow_CO2 = compute_carbon_dioxide_emissions(1, conversion_efficiency, x, MW, LHV, CR)
    specific_emissions = mass_flow_CO2
    return specific_emissions * 1000 * (1000 * 3600)   # Units: g/kWh


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


float_scalar  = nb.typeof(0.)
float_array   = nb.typeof(np.asarray([0.]))
string_scalar = nb.typeof("")
string_array  = nb.typeof(np.asarray([""], dtype='U3'))

spec = [('name',       string_scalar),
        ('components', string_array),
        ('x',          float_array),
        ('y',          float_array),
        ('MW',         float_array),
        ('LHV',        float_array),
        ('CR',         float_array)]

# @jitclass(spec)
class Fluid:

    def __init__(self, name, components, fractions, fraction_type, MW, LHV, CR):

        # Store input variables as instance variables
        self.name       = name
        self.components = components
        self.MW         = MW.astype(np.float64)
        self.LHV        = LHV.astype(np.float64)
        self.CR         = CR.astype(np.float64)

        # Define the fluid composition
        if fraction_type == "molar":
            self.x = fractions
            self.y = convert_molar_to_mass_fraction(self.x, MW)
        elif fraction_type == "mass":
            self.y = fractions
            self.x = convert_mass_to_molar_fraction(self.y, MW)
        else:
            raise Exception("Invalid option for 'fraction_type'. Valid options: 'molar', 'mass'")

        print("hello there!")


        # , fraction_type = "molar"
        # self.fraction_type = fraction_type


        # self.IN_dict    = IN_dict
        # self.name       = IN_dict["name"]         # Fluid name
        # self.MW         = IN_dict["MW"]           # Molecular mass of each component
        # self.LHV        = IN_dict["LHV"]          # Lower heating value of each component
        # self.CR         = IN_dict["CR"]           # Stoichiometric carbon ratio of each component
        # self.components = IN_dict["components"]   # Fluid components

        # # Define the fluid composition
        # if "x" in IN_dict:
        #     self.x = IN_dict["x"]
        #     self.y = convert_molar_to_mass_fraction(self.x, self.MW)
        # elif "y" in IN_dict:
        #     self.y = IN_dict["y"]
        #     self.x = convert_mass_to_molar_fraction(self.y, self.MW)
        # else:
        #     raise Exception("Specify the molar or mass composition of the fluid")
        #
        # # Remove duplicate entries
        # self.cleanup_duplicates()
        #
        # # Compute the average molecular mass
        # self.MW_mean = np.sum(self.x * self.MW)
        #
        # # Compute the average lower heating value
        # self.LHV_mean = np.sum(self.y * self.LHV)


    # def __str__(self):
    #     output = "\n"
    #     output +=  "{} {}\n".format("Fluid name:", self.name)
    #     output += ("{:30}" + "{:>8}"*len(self.components) + "{:>16}" + "\n").format("Components:", *self.components, "Average value")
    #     output += ("{:30}" + "{:>8.2f}" * len(self.components) + "{:>16.2f}" + "\n").format("Molecular mass (g/mol):", *self.MW*1e3, self.MW_mean*1e3)
    #     output += ("{:30}" + "{:>8.2f}" * len(self.components) + "{:>16.2f}" + "\n").format("Heating value (MJ/kg):", *self.LHV/1e6, self.LHV_mean/1e6)
    #     output += ("{:30}" + "{:>8.2f}" * len(self.components) + "\n").format("Carbon ratio (-):", *self.CR)
    #     output += ("{:30}" + "{:>8.4f}" * len(self.components) + "\n").format("Molar composition (-):", *self.x)
    #     output += ("{:30}" + "{:>8.4f}" * len(self.components) + "\n").format("Mass composition (-):",  *self.y)
    #     return output

    # def cleanup_duplicates(self):
    #
    #     # Iterate over the component names
    #     for name in self.components:
    #
    #         # Find the entries that are duplicated
    #         repeated_indices = [i for i, name_bis in enumerate(self.components) if name == name_bis]
    #
    #         # Eliminate repeated information
    #         self.components = np.delete(self.components, repeated_indices[1::])
    #         self.MW   = np.delete(self.MW, repeated_indices[1::])
    #         self.LHV  = np.delete(self.LHV, repeated_indices[1::])
    #
    #         # Adjust the molar fraction of repeated entries
    #         self.x[repeated_indices[0]] = self.x[repeated_indices[0]] + np.sum(self.x[repeated_indices[1::]])
    #         self.x = np.delete(self.x, repeated_indices[1::])
    #
    #         # Adjust the mass fraction of repeated entries
    #         self.y[repeated_indices[0]] = self.y[repeated_indices[0]] + np.sum(self.y[repeated_indices[1::]])
    #         self.y = np.delete(self.y, repeated_indices[1::])






# class Fluid:
#
#     def __init__(self, IN_dict):
#
#         # Store input variables as instance variables
#         self.IN_dict    = IN_dict
#         self.name       = IN_dict["name"]         # Fluid name
#         self.MW         = IN_dict["MW"]           # Molecular mass of each component
#         self.LHV        = IN_dict["LHV"]          # Lower heating value of each component
#         self.CR         = IN_dict["CR"]           # Stoichiometric carbon ratio of each component
#         self.components = IN_dict["components"]   # Fluid components
#
#         # Define the fluid composition
#         if "x" in IN_dict:
#             self.x = IN_dict["x"]
#             self.y = convert_molar_to_mass_fraction(self.x, self.MW)
#         elif "y" in IN_dict:
#             self.y = IN_dict["y"]
#             self.x = convert_mass_to_molar_fraction(self.y, self.MW)
#         else:
#             raise Exception("Specify the molar or mass composition of the fluid")
#
#         # Remove duplicate entries
#         self.cleanup_duplicates()
#
#         # Compute the average molecular mass
#         self.MW_mean = np.sum(self.x * self.MW)
#
#         # Compute the average lower heating value
#         self.LHV_mean = np.sum(self.y * self.LHV)
#
#
#     def __str__(self):
#         output = "\n"
#         output +=  "{} {}\n".format("Fluid name:", self.name)
#         output += ("{:30}" + "{:>8}"*len(self.components) + "{:>16}" + "\n").format("Components:", *self.components, "Average value")
#         output += ("{:30}" + "{:>8.2f}" * len(self.components) + "{:>16.2f}" + "\n").format("Molecular mass (g/mol):", *self.MW*1e3, self.MW_mean*1e3)
#         output += ("{:30}" + "{:>8.2f}" * len(self.components) + "{:>16.2f}" + "\n").format("Heating value (MJ/kg):", *self.LHV/1e6, self.LHV_mean/1e6)
#         output += ("{:30}" + "{:>8.2f}" * len(self.components) + "\n").format("Carbon ratio (-):", *self.CR)
#         output += ("{:30}" + "{:>8.4f}" * len(self.components) + "\n").format("Molar composition (-):", *self.x)
#         output += ("{:30}" + "{:>8.4f}" * len(self.components) + "\n").format("Mass composition (-):",  *self.y)
#         return output
#
#     def cleanup_duplicates(self):
#
#         # Iterate over the component names
#         for name in self.components:
#
#             # Find the entries that are duplicated
#             repeated_indices = [i for i, name_bis in enumerate(self.components) if name == name_bis]
#
#             # Eliminate repeated information
#             self.components = np.delete(self.components, repeated_indices[1::])
#             self.MW   = np.delete(self.MW, repeated_indices[1::])
#             self.LHV  = np.delete(self.LHV, repeated_indices[1::])
#
#             # Adjust the molar fraction of repeated entries
#             self.x[repeated_indices[0]] = self.x[repeated_indices[0]] + np.sum(self.x[repeated_indices[1::]])
#             self.x = np.delete(self.x, repeated_indices[1::])
#
#             # Adjust the mass fraction of repeated entries
#             self.y[repeated_indices[0]] = self.y[repeated_indices[0]] + np.sum(self.y[repeated_indices[1::]])
#             self.y = np.delete(self.y, repeated_indices[1::])


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##

# # Define the properties and composition of natural gas
# natural_gas = create_fluid(
#                 name="Natural gas",
#                 components=np.asarray(["CH4", "C2H6", "C3H8", "C4H10", "C5H12", "C6H14", "N2", "CO2"]),
#                 MW=np.asarray([16.04, 30.07, 44.10, 58.12, 72.15, 86.18, 28.01, 44.01]) / 1e3,
#                 LHV=np.asarray([55.50, 51.90, 50.35, 49.50, 48.60, 47.7, 0.00, 0.00]) * 1e6,
#                 CR=np.asarray([1., 2., 3., 4., 5., 6., 0., 1.]),
#                 fraction=np.asarray([0.72890, 0.13601, 0.08272, 0.02750, 0.00667, 0.00091, 0.01562, 0.00167]),
#                 fraction_type="molar")
#
# # Define properties and composition of methane
# methane = create_fluid(
#                 name="Methane",
#                 components=np.asarray(["CH4"]),
#                 MW=np.asarray([16.04]) / 1e3,
#                 LHV=np.asarray([55.50]) * 1e6,
#                 CR=np.asarray([1]),
#                 fraction=np.asarray([1]),
#                 fraction_type="molar")
#
# # Define properties and composition of hydrogen
# hydrogen = create_fluid(
#                 name="Hydrogen",
#                 components=np.asarray(["H2"]),
#                 MW=np.asarray([2.02]) / 1e3,
#                 LHV=np.asarray([119.96]) * 1e6,
#                 CR=np.asarray([0]),
#                 fraction=np.asarray([1]),
#                 fraction_type="molar")
#
# # Define the properties and composition of carbon monoxide
# carbon_monoxide = create_fluid(
#                 name="Carbon monoxide",
#                 components=np.asarray(["CO"]),
#                 MW=np.asarray([28.01]) / 1e3,
#                 LHV=np.asarray([10.11]) * 1e6,
#                 CR=np.asarray([1]),
#                 fraction=np.asarray([1]),
#                 fraction_type="molar")
#
# # Define the properties and composition of carbon dioxide
# carbon_dioxide = create_fluid(
#                 name="Carbon dioxide",
#                 components=np.asarray(["CO2"]),
#                 MW=np.asarray([44.01]) / 1e3,
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([1]),
#                 fraction=np.asarray([1]),
#                 fraction_type="molar")
#
# # Define the properties and composition of nitrogen
# nitrogen = create_fluid(
#                 name="Nitrogen",
#                 components=np.asarray(["N2"]),
#                 MW=np.asarray([28.01]) / 1e3,
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([0]),
#                 fraction=np.asarray([1]),
#                 fraction_type="molar")
#
# # Define the properties and composition of oxygen
# oxygen = create_fluid(
#                 name="Oxygen",
#                 components=np.asarray(["O2"]),
#                 MW=np.asarray([32.00]) / 1e3,
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([0]),
#                 fraction=np.asarray([1]),
#                 fraction_type="molar")
#
# # Define the properties and composition of water
# water = create_fluid(
#                 name="Water",
#                 components=np.asarray(["H2O"]),
#                 MW=np.asarray([18.02]) / 1e3,
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([0]),
#                 fraction=np.asarray([1]),
#                 fraction_type="molar")
#
# # Define the properties and composition of argon
# argon = create_fluid(
#                 name="Argon",
#                 components=np.asarray(["Ar"]),
#                 MW=np.asarray([39.95]) / 1e3,
#                 LHV=np.asarray([0.00]) * 1e6,
#                 CR=np.asarray([0]),
#                 fraction=np.asarray([1]),
#                 fraction_type="molar")
#
# # Define the properties and composition of air
# air = create_fluid(
#                 name="Air",
#                 components=np.asarray(["N2", "O2", "CO2", "H2O", "Ar"]),
#                 MW=np.asarray([28.01, 32.00, 44.01, 18.02, 39.95]) / 1e3,
#                 LHV=np.asarray([0.00, 0.00, 0.00, 0.00, 0.00]) * 1e6,
#                 CR=np.asarray([0, 0, 1, 0, 0]),
#                 fraction=np.asarray([0.7739, 0.2076, 0.0003, 0.0089, 0.0093]),
#                 fraction_type="molar")
