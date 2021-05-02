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
import os
import pdb
import time
import numpy as np
import matplotlib.pyplot as plt

import hes_off_core


# # Evaluate the process model (functional with Numba)
# IN = hes_off.read_configuration_file("edward_grieg_GT.cfg")
# EdwardGrieg = hes_off.IntegratedModel(IN)
# with hes_off.Timer():
#     EdwardGrieg.evaluate_process_model()
# print("Lifetime CO2 emissions {:0.2f} Mton".format(EdwardGrieg.CO2_emissions_total / 1e9))


# IN = hes_off.read_configuration_file("edward_grieg_GT_WT.cfg")
# EdwardGrieg = hes_off.IntegratedModel(IN)
# with hes_off.Timer():
#     EdwardGrieg.evaluate_process_model()
# print("Lifetime CO2 emissions {:0.2f} Mton".format(EdwardGrieg.CO2_emissions_total / 1e9))


# EdwardGrieg.plot_hydrogen_level()
# EdwardGrieg.plot_hydrogen_balance()

IN = hes_off_core.read_configuration_file("edward_grieg.cfg")
EdwardGrieg = hes_off_core.IntegratedModel(IN)
with hes_off_core.Timer():
    EdwardGrieg.evaluate_process_model()
# print("Lifetime CO2 emissions {:0.2f} Mton".format(EdwardGrieg.CO2_emissions / 1e9))
print(EdwardGrieg.CO2_emissions/1e9)


a =EdwardGrieg.plot_hydrogen_level()
EdwardGrieg.plot_hydrogen_balance()
EdwardGrieg.plot_power_balance()
EdwardGrieg.plot_power_deficit()
EdwardGrieg.plot_carbon_dioxide_emissions()
# a.savefig("test.png")
plt.show()




# EdwardGrieg.plot_sensitivity_analysis(variable="WT_RATED_POWER",
#                                        low_value=EdwardGrieg.WT_RATED_POWER*0.00,
#                                        high_value=EdwardGrieg.WT_RATED_POWER*2.00)
#
# EdwardGrieg.plot_sensitivity_analysis(variable="H2_CAPACITY",
#                                        low_value=EdwardGrieg.H2_CAPACITY*0.50,
#                                        high_value=EdwardGrieg.H2_CAPACITY*2.00)
# #
# EdwardGrieg.plot_sensitivity_analysis(variable="GT_MAX_H2",
#                                        low_value=0.0,
#                                        high_value=1.00)

# # EdwardGrieg.plot_sensitivity_analysis(variable="EL_RATED_POWER",
# #                                        low_value=EdwardGrieg.EL_RATED_POWER*0.50,
# #                                        high_value=EdwardGrieg.EL_RATED_POWER*2.00)
# #
# # EdwardGrieg.plot_sensitivity_analysis(variable="FC_RATED_POWER",
# #                                        low_value=EdwardGrieg.FC_RATED_POWER*0.50,
# #                                        high_value=EdwardGrieg.FC_RATED_POWER*2.00)
#


# # Plot the process results
# EdwardGrieg.plot_hydrogen_level()
# EdwardGrieg.plot_hydrogen_balance()
# EdwardGrieg.plot_power_deficit()
# EdwardGrieg.plot_power_balance()
# EdwardGrieg.plot_carbon_dioxide_emissions()
# EdwardGrieg.plot_flag()



# #
#
#
# plt.show()
# #
# #
# #
#
#
#
#
# #     # for i in range(30):
# #     #     with hes_off_object_oriented.Timer():
# #     #
# #     #         FC_power = FC_RATED_POWER * 0.9
# #     #         # FC_power = np.asarray([FC_RATED_POWER * 0.9, FC_RATED_POWER * 0.9])
# #     #         # compute_FC_efficiency(model='POWERCELL_S3', efficiency_coefficients=np.asarray([0]), rated_power=FC_RATED_POWER, power_output=FC_RATED_POWER * 0.9)
# #     #         # compute_FC_specific_performance(model='POWERCELL_S3', efficiency_coefficients=0.0, rated_power=FC_RATED_POWER, power_output=FC_RATED_POWER * 0.9)
# #     #         # compute_FC_hydrogen_consumption(model='POWERCELL_S3', efficiency_coefficients=FC_EFFICIENCY, rated_power=FC_RATED_POWER, power_output=FC_power)
# #     #
# #     #
# #     #         # compute_EL_efficiency(model='NEL_HYDROGEN', efficiency_coefficients=0.0, rated_power=EL_RATED_POWER, power_input=EL_RATED_POWER * 0.9)
# #     #         # compute_EL_specific_performance(model='NEL_HYDROGEN', efficiency_coefficients=0.0, rated_power=EL_RATED_POWER, power_input=EL_RATED_POWER * 0.9)
# #     #         # compute_EL_hydrogen_production(model='NEL_HYDROGEN', efficiency_coefficients=EL_EFFICIENCY, rated_power=EL_RATED_POWER, power_input=EL_RATED_POWER * 0.9)
# #     #         #
# #     #         # compute_GT_efficiency(model='LM2500+G4', number_of_units=1, power_output=30e6)
# #     #         # compute_GT_power_from_heat(model='LM2500+G4', number_of_units=1, heat_output=20e6)
# #     #         # compute_GT_power_from_heat(model='LM6000-PF', number_of_units=1, heat_output=20e6)
# #     #
# #     #         # compute_WT_power_output(model="HYWIND", rated_power=WT_RATED_POWER, hub_height=WT_HUB_HEIGHT, ref_height=WT_HUB_HEIGHT, wind_speed=wind_speed)
# #     #         # compute_WT_capacity_factor(model="HYWIND", rated_power=30e6, hub_height=90, ref_height=14, wind_speed=wind_speed, time_series=wind_time)
# #     #
# #     #         evaluate_process_model(HEAT_DEMAND, POWER_DEMAND,
# #     #                                GT_MODEL, GT_UNITS, GT_MAX_H2,
# #     #                                WT_MODEL, WT_RATED_POWER, WT_REF_HEIGHT, WT_HUB_HEIGHT,
# #     #                                EL_MODEL, EL_RATED_POWER, EL_EFFICIENCY,
# #     #                                FC_MODEL, FC_RATED_POWER, FC_EFFICIENCY,
# #     #                                H2_CAPACITY, H2_INITIAL_FRACTION, H2_RECHARGE_THRESHOLD, H2_COFIRE_THRESHOLD,
# #     #                                wind_speed)
# #
