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
import hes_off_object_oriented
import numpy as np
import matplotlib.pyplot as plt

# Defined wind turbine tests
def test_load_wind_data():

    # Load wind datafile
    filename = 'SLEIPNERWIND'
    windData = hes_off_object_oriented.process_models.read_wind_data(filename=filename)

    # Check that the wind datafile was loaded correctly
    assert np.abs(np.min(windData["speed"]  - 0.000000000000)) < 1e-6
    assert np.abs(np.mean(windData["speed"] - 7.815201701644)) < 1e-6
    assert np.abs(np.max(windData["speed"]  - 29.29416666572)) < 1e-6


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


def test_WT_initialization():

    # Empty initialization
    testWT = hes_off_object_oriented.process_models.WT()
    print(testWT)

    # HYWIND turbine
    testWT = hes_off_object_oriented.process_models.WT(model='HYWIND', rated_power=30e6, wind_speed_factor=1.25)
    print(testWT)

    # NREL turbine
    testWT = hes_off_object_oriented.process_models.WT(model='NREL', rated_power=30e6, wind_speed_factor=(90, 14))
    print(testWT)

def test_WT_power_computation():

    # HYWIND turbines at 20 m/s
    testWT = hes_off_object_oriented.process_models.WT(model='HYWIND', wind_speed_factor=1)
    power = testWT.compute_power_output(wind_speed=20)
    assert testWT._rated_power == power

    # NREL turbines at 20 m/s
    testWT = hes_off_object_oriented.process_models.WT(model='NREL', wind_speed_factor=1)
    power = testWT.compute_power_output(wind_speed=20)
    assert testWT._rated_power == power

def test_WT_capacity_factor_computation_1():

    # One year with 20 m/s wind speed
    time = np.arange(1, 8760+2)
    speed = 20 + 0*time

    # HYWIND turbine
    testWT = hes_off_object_oriented.process_models.WT(model='HYWIND', wind_speed_factor=1)
    capacity_factor = testWT.compute_capacity_factor(wind_speed=speed, time_series=time)
    assert capacity_factor == 1.0

    # NREL turbine
    testWT = hes_off_object_oriented.process_models.WT(model='NREL', wind_speed_factor=1)
    capacity_factor = testWT.compute_capacity_factor(wind_speed=speed, time_series=time)
    assert capacity_factor == 1.0

def test_WT_capacity_factor_computation_2():

    # Load wind datafile
    filename = 'SLEIPNERWIND'
    windData = hes_off_object_oriented.process_models.read_wind_data(filename=filename)

    # HYWIND turbine
    testWT = hes_off_object_oriented.process_models.WT(model='HYWIND', wind_speed_factor=(90, 14))
    capacity_factor = testWT.compute_capacity_factor(wind_speed=windData['speed'], time_series=windData["time"])
    assert np.abs(capacity_factor - 0.5528896376081194) < 1e-9      # Value computed in advance

    # NREL turbine
    testWT = hes_off_object_oriented.process_models.WT(model='NREL', wind_speed_factor=(90, 14))
    capacity_factor = testWT.compute_capacity_factor(wind_speed=windData['speed'], time_series=windData["time"])
    assert np.abs(capacity_factor - 0.6565242884811453) < 1e-9      # Value computed in advance


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


def test_GT_initialization():

    # Empty initialization
    testGT = hes_off_object_oriented.process_models.GT()
    print(testGT)

    # LM2500+G4 initialization
    testGT = hes_off_object_oriented.process_models.GT(model='LM2500+G4')
    print(testGT)

    # Update the number of gas turbines after initialization
    testGT.number_of_units = 3
    print(testGT)

    # LM6000-PF initialization
    testGT = hes_off_object_oriented.process_models.GT(model='LM6000-PF', number_of_units=2)
    print(testGT)

# def test_GT_power_output():
#
#     # Test the LM2500+G4 model
#     min_load, max_load = 0.00, 1.00
#     mid_load = (min_load + max_load)/2
#     testGT = hes_off_object_oriented.process_models.GT(model='LM2500+G4', min_load=min_load, max_load=max_load,number_of_units=1)
#     assert np.abs(testGT.compute_power_output(-1.0*mid_load*testGT.rated_power) - 0.00) < 1e-9
#     assert np.abs(testGT.compute_power_output(+0.5*min_load*testGT.rated_power) - min_load*testGT.rated_power) < 1e-9
#     assert np.abs(testGT.compute_power_output(+1.0*mid_load*testGT.rated_power) - mid_load*testGT.rated_power) < 1e-9
#     assert np.abs(testGT.compute_power_output(+2.0*max_load*testGT.rated_power) - max_load*testGT.rated_power) < 1e-9
#
#     # Test the LM6000-PF model
#     min_load, max_load = 0.10, 0.90
#     mid_load = (min_load + max_load)/2
#     testGT = hes_off_object_oriented.process_models.GT(model='LM6000-PF', min_load=min_load, max_load=max_load, number_of_units=2)
#     assert np.abs(testGT.compute_power_output(-1.0*mid_load*testGT.rated_power) - 0.00) < 1e-9
#     assert np.abs(testGT.compute_power_output(+0.5*min_load*testGT.rated_power) - min_load*testGT.rated_power) < 1e-9
#     assert np.abs(testGT.compute_power_output(+1.0*mid_load*testGT.rated_power) - mid_load*testGT.rated_power) < 1e-9
#     assert np.abs(testGT.compute_power_output(+2.0*max_load*testGT.rated_power) - max_load*testGT.rated_power) < 1e-9

def test_GT_efficiency():

    # Test the LM2500+G4 model
    testGT = hes_off_object_oriented.process_models.GT(model='LM2500+G4', number_of_units=1)
    assert np.abs(testGT.compute_efficiency(testGT.rated_power) - testGT.design_efficiency) < 1e-9

    # Test the LM6000-PF model
    testGT = hes_off_object_oriented.process_models.GT(model='LM6000-PF', number_of_units=2)
    assert np.abs(testGT.compute_efficiency(testGT.rated_power) - testGT.design_efficiency) < 1e-9

def test_GT_emissions_1():

    # Define the fluids
    CH4 = hes_off_object_oriented.combustion_models.methane
    CO2 = hes_off_object_oriented.combustion_models.carbon_dioxide

    # Compute the CO2 emissions of the LM2500+G4 model at full load
    testGT = hes_off_object_oriented.process_models.GT(model='LM2500+G4', number_of_units=1)
    CO2_flow_rate = testGT.compute_carbon_dioxide_emissions(power_output=testGT.rated_power, fuel=CH4)
    CO2_flow_rate_bis = (CO2.MW / CH4.MW) * (1 / CH4.LHV) * (testGT.rated_power / testGT.design_efficiency)
    CO2_specific = testGT.compute_specific_carbon_dioxide_emissions(power_output=testGT.rated_power, fuel=CH4)
    CO2_specific_bis = CO2_flow_rate_bis / testGT.rated_power * 1e3 * 3.6e6
    assert np.abs(CO2_flow_rate - CO2_flow_rate_bis) < 1e-9
    assert np.abs(CO2_specific - CO2_specific_bis) < 1e-9

    # Compute the CO2 emissions of the LM6000-PF model at full load
    testGT = hes_off_object_oriented.process_models.GT(model='LM6000-PF', number_of_units=1)
    CO2_flow_rate = testGT.compute_carbon_dioxide_emissions(power_output=testGT.rated_power, fuel=CH4)
    CO2_flow_rate_bis = (CO2.MW / CH4.MW) * (1 / CH4.LHV) * (testGT.rated_power / testGT.design_efficiency)
    CO2_specific = testGT.compute_specific_carbon_dioxide_emissions(power_output=testGT.rated_power, fuel=CH4)
    CO2_specific_bis = CO2_flow_rate_bis / testGT.rated_power * 1e3 * 3.6e6
    assert np.abs(CO2_flow_rate - CO2_flow_rate_bis) < 1e-9
    assert np.abs(CO2_specific - CO2_specific_bis) < 1e-9

def test_GT_emissions_2():

    # Define the fluids
    H2 = hes_off_object_oriented.combustion_models.hydrogen

    # Compute the CO2 emissions of the LM2500+G4 model at full load
    testGT = hes_off_object_oriented.process_models.GT(model='LM2500+G4', number_of_units=1)
    CO2_flow_rate = testGT.compute_carbon_dioxide_emissions(power_output=testGT.rated_power, fuel=H2)
    CO2_specific = testGT.compute_specific_carbon_dioxide_emissions(power_output=testGT.rated_power, fuel=H2)
    assert np.abs(CO2_flow_rate - 0.0) < 1e-9
    assert np.abs(CO2_specific - 0.0) < 1e-9

    # Compute the CO2 emissions of the LM6000-PF model at full load
    testGT = hes_off_object_oriented.process_models.GT(model='LM6000-PF', number_of_units=1)
    CO2_flow_rate = testGT.compute_carbon_dioxide_emissions(power_output=testGT.rated_power, fuel=H2)
    CO2_specific = testGT.compute_specific_carbon_dioxide_emissions(power_output=testGT.rated_power, fuel=H2)
    assert np.abs(CO2_flow_rate - 0.0) < 1e-9
    assert np.abs(CO2_specific - 0.0) < 1e-9

def test_GT_heat_demand():

    # Check that the heat computations are consistent
    testGT = hes_off_object_oriented.process_models.GT(model='LM2500+G4', number_of_units=1)
    heat_demand = testGT.rated_power
    power_output = testGT.compute_minimum_power_output(heat_demand)
    heat_output = testGT.compute_heat_from_power(power_output)
    assert np.abs(heat_demand - heat_output) < 1e-3


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


def test_EL_initialization():

    # Empty initialization
    testEL = hes_off_object_oriented.process_models.EL()
    print(testEL)

    # Polynomial efficiency model (constant efficiency)
    testEL = hes_off_object_oriented.process_models.EL(model='POLYNOMIAL_EFFICIENCY', efficiency_coefficients=0.8, rated_power=5e6)
    print(testEL)

    # Polynomial efficiency model (NEL_HYDROGEN)
    coeff_NEL_HYDROGEN = np.asarray([6.82299e-03, -8.29304e-09, 1.59861e-14, -1.41175e-20]) / 1e6 * hes_off_object_oriented.hydrogen.LHV
    testEL = hes_off_object_oriented.process_models.EL(model='POLYNOMIAL_EFFICIENCY', efficiency_coefficients=coeff_NEL_HYDROGEN, rated_power=5e6)
    testEL.rated_power = 10e6
    testEL.unit_power = 1e6
    print(testEL)

def test_EL_hydrogen_production_computation_1():

    # NEL_HYDROGEN stack
    testEL = hes_off_object_oriented.process_models.EL(model='NEL_HYDROGEN', rated_power=2e6)
    SHP = testEL.compute_specific_hydrogen_production(1.0e6)*1e6
    EFF = testEL.compute_hydrogen_conversion_efficiency(1.0e6)
    assert np.abs(SHP - 0.0055282753125) < 1e-9      # Value computed in advance
    assert np.abs(EFF - 0.6631719064875) < 1e-9      # Value computed in advance
    assert testEL.number_of_units == 4

def test_EL_hydrogen_production_computation_2():

    # Constant efficiency stack
    eff = 0.75
    LHV = hes_off_object_oriented.combustion_models.hydrogen.LHV
    testEL = hes_off_object_oriented.process_models.EL(model='POLYNOMIAL_EFFICIENCY', efficiency_coefficients=eff, rated_power=1e6)
    SHP = testEL.compute_specific_hydrogen_production(1.0e6)
    EFF = testEL.compute_hydrogen_conversion_efficiency(1.0e6)
    assert np.abs(SHP - eff/LHV) < 1e-9      # Value computed in advance
    assert np.abs(EFF - eff)     < 1e-9      # Value computed in advance


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


def test_FC_initialization():

    # Empty initialization
    testFC = hes_off_object_oriented.process_models.FC()
    print(testFC)

    # Polynomial efficiency model
    testFC = hes_off_object_oriented.process_models.FC(model='POLYNOMIAL_EFFICIENCY', efficiency_coefficients=0.6, rated_power=5e6)
    print(testFC)

    # Polynomial specific energy production model
    testFC = hes_off_object_oriented.process_models.FC(model='POLYNOMIAL_EFFICIENCY', efficiency_coefficients=[0.6], rated_power=5e6)
    testFC.rated_power = 10e6
    testFC.unit_power = 1e6
    print(testFC)

def test_FC_energy_production_computation_1():

    # POWERCELL_S3 stack
    testFC = hes_off_object_oriented.process_models.FC(model='POWERCELL_S3', rated_power=2e6)
    power_output = np.linspace(0, testFC.rated_power, 500)
    SEP = testFC.compute_specific_energy_production(power_output)
    assert np.abs(SEP[0]/1e6      - 0)                  < 1e-6  # Value computed in advance
    assert np.abs(SEP[-1]/1e6     - 56.24397168621687)  < 1e-9  # Value computed in advance
    assert np.abs(np.max(SEP)/1e6 - 64.22908330778608)  < 1e-9  # Value computed in advance
    assert testFC.number_of_units == 2e6/125e3

def test_FC_energy_production_computation_2():

    # Constant efficiency stack
    eff = 0.50
    LHV = hes_off_object_oriented.combustion_models.hydrogen.LHV
    testFC = hes_off_object_oriented.process_models.FC(model='POLYNOMIAL_EFFICIENCY', efficiency_coefficients=eff, rated_power=1e6)
    SEP = testFC.compute_specific_energy_production(1.0e6)
    EFF = testFC.compute_hydrogen_conversion_efficiency(1.0e6)
    H2_flow = testFC.compute_hydrogen_consumption(1.0e6)
    assert np.abs(SEP - eff*LHV) < 1e-9        # Value computed in advance
    assert np.abs(EFF - eff)       < 1e-9      # Value computed in advance
    assert np.abs(H2_flow - 1.0e6/LHV/eff) < 1e-9


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


# Run the tests when this script is executed as main
if __name__ == '__main__':
    test_load_wind_data()
    test_WT_initialization()
    test_WT_power_computation()
    test_WT_capacity_factor_computation_1()
    test_WT_capacity_factor_computation_2()
    test_GT_initialization()
    # test_GT_power_output()
    test_GT_efficiency()
    test_GT_emissions_1()
    test_GT_emissions_2()
    test_GT_heat_demand()
    test_EL_initialization()
    test_EL_hydrogen_production_computation_1()
    test_EL_hydrogen_production_computation_2()
    test_FC_initialization()
    test_FC_energy_production_computation_1()
    test_FC_energy_production_computation_2()



