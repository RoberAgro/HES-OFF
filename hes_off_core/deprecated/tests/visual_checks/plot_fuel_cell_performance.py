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
import matplotlib.pyplot as plt

# Create output folder if it does not exist
if os.path.exists('figures/') is False:
    os.mkdir('figures')

# Define font settings
plt.rc('text', usetex=False)
plt.rcParams['font.family']      = 'serif'      # 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'
plt.rcParams['font.serif']       = 'cmr10'      # 'cmr10', 'palatino', 'times new roman'
plt.rcParams['mathtext.fontset'] = 'cm'         # 'cm' (latex style), 'stix' (times new roman style), 'stixsans'

# Plot NEL Hydrogen performance curves
FC = hes_off_object_oriented.process_models.FC(model='POWERCELL_S3', rated_power=125e3)
fig1, _ = FC.plot_hydrogen_consumption_curve()
fig2, _ = FC.plot_specific_electricity_production_curve()
fig3, _ = FC.plot_hydrogen_conversion_efficiency_curve()
fig1.savefig('figures/PowerCell_performance.pdf', bbox_inches='tight')
fig2.savefig('figures/PowerCell_specific_performance.pdf', bbox_inches='tight')
fig3.savefig('figures/PowerCell_efficiency.pdf', bbox_inches='tight')

# Show the figure
plt.show()






