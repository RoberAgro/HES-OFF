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
import numpy as np


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


def convert_molar_to_mass_fraction(x, MW):
    x, MW = np.asarray(x), np.asarray(MW)
    if np.abs(1.0 - np.sum(x)) > 1e-12:
        raise Exception("The sum of the input molar fractions is not unity")
    y = (x*MW)/np.sum(x*MW)
    return y

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
    MW_CO2 = carbon_dioxide.MW
    mass_flow_CO2 = np.sum(CR * x, axis=0) / np.sum(x * MW * LHV, axis=0) * (power_output / conversion_efficiency) * MW_CO2
    return mass_flow_CO2   # Units: kg/s

def compute_specific_carbon_dioxide_emissions(conversion_efficiency, x, MW, LHV, CR):
    mass_flow_CO2 = compute_carbon_dioxide_emissions(1, conversion_efficiency, x, MW, LHV, CR)
    specific_emissions = mass_flow_CO2
    return specific_emissions * 1000 * (1000 * 3600)   # Units: g/kWh

def create_fluid_mixture(fluids, fractions, fraction_type="molar", name="Fuel mixture"):

    # Check if the composition fractions are consistent
    if np.abs(1.0 - np.sum(fractions)) > 1e-12:
        raise Exception("The sum of the input fractions is not unity")

    # Create fluid mixture
    mixture = {}
    mixture["name"] = name
    mixture["MW"]  = np.concatenate([fluid.MW   for fluid in fluids])
    mixture["LHV"] = np.concatenate([fluid.LHV  for fluid in fluids])
    mixture["CR"]  = np.concatenate([fluid.CR   for fluid in fluids])
    mixture["components"] = np.concatenate([fluid.components for fluid in fluids])
    if fraction_type == "molar":
        mixture["x"] = np.concatenate([fraction*fluid.x for fraction, fluid in zip(fractions, fluids)])
        mixture["y"] = convert_molar_to_mass_fraction(mixture["x"], mixture["MW"])
    elif fraction_type == "mass":
        mixture["y"] = np.concatenate([fraction*fluid.y for fraction, fluid in zip(fractions, fluids)])
        mixture["x"] = convert_mass_to_molar_fraction(mixture["y"], mixture["MW"])
    else:
        raise Exception("Invalid option for 'fraction_type'. Valid options: 'molar', 'mass'")

    return Fluid(mixture)


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


class Fluid:

    def __init__(self, IN_dict):

        # Store input variables as instance variables
        self.IN_dict    = IN_dict
        self.name       = IN_dict["name"]         # Fluid name
        self.MW         = IN_dict["MW"]           # Molecular mass of each component
        self.LHV        = IN_dict["LHV"]          # Lower heating value of each component
        self.CR         = IN_dict["CR"]           # Stoichiometric carbon ratio of each component
        self.components = IN_dict["components"]   # Fluid components

        # Define the fluid composition
        if "x" in IN_dict:
            self.x = IN_dict["x"]
            self.y = convert_molar_to_mass_fraction(self.x, self.MW)
        elif "y" in IN_dict:
            self.y = IN_dict["y"]
            self.x = convert_mass_to_molar_fraction(self.y, self.MW)
        else:
            raise Exception("Specify the molar or mass composition of the fluid")

        # Remove duplicate entries
        self.cleanup_duplicates()

        # Compute the average molecular mass
        self.MW_mean = np.sum(self.x * self.MW)

        # Compute the average lower heating value
        self.LHV_mean = np.sum(self.y * self.LHV)


    def __str__(self):
        output = "\n"
        output +=  "{} {}\n".format("Fluid name:", self.name)
        output += ("{:30}" + "{:>8}"*len(self.components) + "{:>16}" + "\n").format("Components:", *self.components, "Average value")
        output += ("{:30}" + "{:>8.2f}" * len(self.components) + "{:>16.2f}" + "\n").format("Molecular mass (g/mol):", *self.MW*1e3, self.MW_mean*1e3)
        output += ("{:30}" + "{:>8.2f}" * len(self.components) + "{:>16.2f}" + "\n").format("Heating value (MJ/kg):", *self.LHV/1e6, self.LHV_mean/1e6)
        output += ("{:30}" + "{:>8.2f}" * len(self.components) + "\n").format("Carbon ratio (-):", *self.CR)
        output += ("{:30}" + "{:>8.4f}" * len(self.components) + "\n").format("Molar composition (-):", *self.x)
        output += ("{:30}" + "{:>8.4f}" * len(self.components) + "\n").format("Mass composition (-):",  *self.y)
        return output

    def cleanup_duplicates(self):

        # Iterate over the component names
        for name in self.components:

            # Find the entries that are duplicated
            repeated_indices = [i for i, name_bis in enumerate(self.components) if name == name_bis]

            # Eliminate repeated information
            self.components = np.delete(self.components, repeated_indices[1::])
            self.MW   = np.delete(self.MW, repeated_indices[1::])
            self.LHV  = np.delete(self.LHV, repeated_indices[1::])

            # Adjust the molar fraction of repeated entries
            self.x[repeated_indices[0]] = self.x[repeated_indices[0]] + np.sum(self.x[repeated_indices[1::]])
            self.x = np.delete(self.x, repeated_indices[1::])

            # Adjust the mass fraction of repeated entries
            self.y[repeated_indices[0]] = self.y[repeated_indices[0]] + np.sum(self.y[repeated_indices[1::]])
            self.y = np.delete(self.y, repeated_indices[1::])


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


# Define properties and composition of natural gas
fluid_dict = {"name":       "Natural gas",
              "components": np.asarray(["CH4", "C2H6", "C3H8", "C4H10", "C5H12", "C6H14", "N2", "CO2"]),
              "MW":         np.asarray([16.04, 30.07, 44.10, 58.12, 72.15, 86.18, 28.01, 44.01]) / 1e3,
              "LHV":        np.asarray([55.50, 51.90, 50.35, 49.50, 48.60, 47.7, 0.00, 0.00]) * 1e6,
              "CR":         np.asarray([1, 2, 3, 4, 5, 6, 0, 1]),
              "x":          np.asarray([0.72890, 0.13601, 0.08272, 0.02750, 0.00667, 0.00091, 0.01562, 0.00167])}
natural_gas = Fluid(fluid_dict)

# Define properties and composition of methane
fluid_dict = {"name":       "Methane",
              "components": np.asarray(["CH4"]),
              "MW":         np.asarray([16.04]) / 1e3,
              "LHV":        np.asarray([55.50]) * 1e6,
              "CR":         np.asarray([1]),
              "x":          np.asarray([1])}
methane = Fluid(fluid_dict)

# Define the properties and composition of hydrogen
fluid_dict = {"name":      "Hydrogen",
              "components": np.asarray(["H2"]),
              "MW":         np.asarray([2.02]) / 1e3,
              "LHV":        np.asarray([119.96]) * 1e6,
              "CR":         np.asarray([0]),
              "x":          np.asarray([1])}
hydrogen = Fluid(fluid_dict)

# Define the properties and composition of carbon monoxide
fluid_dict = {"name":       "Carbon monoxide",
              "components": np.asarray(["CO"]),
              "MW":         np.asarray([28.01]) / 1e3,
              "LHV":        np.asarray([10.11]) * 1e6,
              "CR":         np.asarray([1]),
              "y":          np.asarray([1])}
carbon_monoxide = Fluid(fluid_dict)


# Define the properties and composition of carbon dioxide
fluid_dict = {"name":       "Carbon dioxide",
              "components": np.asarray(["CO2"]),
              "MW":         np.asarray([44.01]) / 1e3,
              "LHV":        np.asarray([0.00]) * 1e6,
              "CR":         np.asarray([1]),
              "y":          np.asarray([1])}
carbon_dioxide = Fluid(fluid_dict)

# Define the properties and composition of air
fluid_dict = {"name":       "Air",
              "components": np.asarray(["N2", "O2", "CO2", "H2O", "Ar"]),
              "MW":         np.asarray([28.01, 32.00, 44.01, 18.02, 39.95]) / 1e3,
              "LHV":        np.asarray([0.00, 0.00, 0.00, 0.00, 0.00]) * 1e6,
              "CR":         np.asarray([0, 0, 1, 0, 0]),
              "x":          np.asarray([0.7739, 0.2076, 0.0003, 0.0089, 0.0093])}
air = Fluid(fluid_dict)

# Define the properties and composition of nitrogen
fluid_dict = {"name":       "Nitrogen",
              "components": np.asarray(["N2"]),
              "MW":         np.asarray([28.01]) / 1e3,
              "LHV":        np.asarray([0.00]) * 1e6,
              "CR":         np.asarray([0]),
              "x":          np.asarray([1])}
nitrogen = Fluid(fluid_dict)

# Define the properties and composition of oxygen
fluid_dict = {"name":       "Oxygen",
              "components": np.asarray(["O2"]),
              "MW":         np.asarray([32.00]) / 1e3,
              "LHV":        np.asarray([0.00]) * 1e6,
              "CR":         np.asarray([0]),
              "x":          np.asarray([1])}
oxygen = Fluid(fluid_dict)

# Define the properties and composition of water
fluid_dict = {"name":       "Water",
              "components": np.asarray(["H20"]),
              "MW":         np.asarray([18.02]) / 1e3,
              "LHV":        np.asarray([0.00]) * 1e6,
              "CR":         np.asarray([0]),
              "x":          np.asarray([1])}
water = Fluid(fluid_dict)

# Define the properties and composition of argon
fluid_dict = {"name":       "Argon",
              "components": np.asarray(["Ar"]),
              "MW":         np.asarray([39.95]) / 1e3,
              "LHV":        np.asarray([0.00]) * 1e6,
              "CR":         np.asarray([0]),
              "x":          np.asarray([1])}
argon = Fluid(fluid_dict)

# Define the properties and composition of coal
fluid_dict = {"name":       "Coal",
              "components": np.asarray(["Coal"]),
              "MW":         np.asarray([12.00]) / 1e3,
              "LHV":        np.asarray([32.80]) * 1e6,
              "CR":         np.asarray([1]),
              "x":          np.asarray([1.00])}
coal = Fluid(fluid_dict)

