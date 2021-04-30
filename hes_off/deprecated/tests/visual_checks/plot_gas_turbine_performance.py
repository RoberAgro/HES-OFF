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
import hes_off_object_oriented
import numpy as np
import matplotlib.pyplot as plt

# Create output folder if it does not exist
if os.path.exists('figures/') is False:
    os.mkdir('figures')

# Define font settings
plt.rc('text', usetex=False)
plt.rcParams['font.family']      = 'serif'      # 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'
plt.rcParams['font.serif']       = 'cmr10'      # 'cmr10', 'palatino', 'times new roman'
plt.rcParams['mathtext.fontset'] = 'cm'         # 'cm' (latex style), 'stix' (times new roman style), 'stixsans'

# Plot the efficiency and specific CO2 emissions of the LM2500+G4 gas turbine
GT = hes_off_object_oriented.process_models.GT(model="LM2500+G4", number_of_units=1)
fig1, _ = GT.plot_efficiency_curve()
fig2, _ = GT.plot_specific_carbon_dioxide_emissions_curve()
fig3, _ = GT.plot_heat_and_power_curve()
fig1.savefig('figures/LM2500+G4_efficiency_curve.pdf', bbox_inches='tight')
fig2.savefig('figures/LM2500+G4_CO2_emissions_curve.pdf', bbox_inches='tight')

# Plot the efficiency and specific CO2 emissions of the LM6000-PF gas turbine
GT = hes_off_object_oriented.process_models.GT(model="LM6000-PF", number_of_units=1)
fig4, _ = GT.plot_efficiency_curve()
fig5, _ = GT.plot_specific_carbon_dioxide_emissions_curve()
fig6, _ = GT.plot_heat_and_power_curve()
fig4.savefig('figures/LM6000-PF_efficiency_curve.pdf', bbox_inches='tight')
fig5.savefig('figures/LM6000-PF_CO2_emissions_curve.pdf', bbox_inches='tight')

# Show the figures
plt.show()
