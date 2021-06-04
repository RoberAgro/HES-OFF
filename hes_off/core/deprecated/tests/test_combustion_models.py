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
import matplotlib.pyplot as plt
import hes_off_object_oriented.combustion_models as cm

## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


def test_fluid_initialization():

    # Initialize fluid from user-defined dictionary
    fluid_dict = {"name": "Air",
                  "components": np.asarray(["N2", "O2", "CO2", "H2O", "Ar"]),
                  "MW":         np.asarray([28.01, 32.00, 44.01, 18.02, 39.95]) / 1e3,
                  "LHV":        np.asarray([00.00, 00.00, 00.00, 00.00, 00.00]) * 1e6,
                  "CR":         np.asarray([0, 0, 1, 0, 0]),
                  "x":          np.asarray([0.7739,0.2076,0.0003,0.0089,0.0093])}

    air = cm.Fluid(fluid_dict)
    print(air)

    # Initialize fluid from the input dictionary of an existing fluid
    fluid_dict = cm.natural_gas.IN
    natural_gas = cm.Fluid(fluid_dict)
    print(natural_gas)


def test_mixture_creation():

    # Create an air mixture using a molar basis and check its molecular mass
    components = [cm.nitrogen, cm.oxygen, cm.carbon_dioxide, cm.water, cm.argon]
    molar_composition = [0.7739, 0.2076, 0.0003, 0.0089, 0.0093]
    air1 = cm.create_fluid_mixture(components, molar_composition, fraction_type="molar")
    assert np.abs(air1.MW_mean - 0.028865255) < 1e-9

    # Create an air mixture using a mass basis and check its molecular mass
    air2 = cm.create_fluid_mixture(components, air1.y, fraction_type="mass")
    assert np.abs(air2.MW_mean - 0.028865255) < 1e-9


def test_molar_to_mass_conversion():

    # Simple numerical test that can be solved "by hand"
    x  = [0.5, 0.5]
    MW = [10, 20]
    y  = cm.convert_molar_to_mass_fraction(x, MW)
    assert np.abs(y[0] - x[0]*10/15) < 1e-9
    assert np.abs(y[1] - x[1]*20/15) < 1e-9

def test_mass_to_molar_conversion():

    # Simple numerical test that can be solved "by hand"
    y  = [0.5, 0.5]
    MW = [10, 20]
    x  = cm.convert_mass_to_molar_fraction(y, MW)
    assert np.abs(x[0] - y[0]*(2*20)/(10+20)) < 1e-9
    assert np.abs(x[1] - y[1]*(2*10)/(10+20)) < 1e-9

def test_carbon_dioxide_emissions_computation():

    # Compute the specific CO2 emissions of natural gas combustion
    efficiency = 1.00
    fuel = cm.natural_gas
    specific_emissions = cm.compute_specific_carbon_dioxide_emissions(efficiency, fuel.x, fuel.MW, fuel.LHV, fuel.CR)
    assert np.abs(specific_emissions - 192.6400256508272) < 1e-9   # Value computed in advance

    # Compute the specific CO2 emissions of hydrogen combustion
    efficiency = 1.00
    fuel = cm.hydrogen
    specific_emissions = cm.compute_specific_carbon_dioxide_emissions(efficiency, fuel.x, fuel.MW, fuel.LHV, fuel.CR)
    assert np.abs(specific_emissions - 0.00) < 1e-9   # Should be zero


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


# Run the tests when this script is executed as main
if __name__ == '__main__':
    test_fluid_initialization()
    test_mixture_creation()
    test_molar_to_mass_conversion()
    test_mass_to_molar_conversion()
    test_carbon_dioxide_emissions_computation()