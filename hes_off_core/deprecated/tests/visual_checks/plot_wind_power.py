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

# Plot HYWIND power curve
WT = hes_off_object_oriented.process_models.WT(model='HYWIND')
fig, _ = WT.plot_power_curve()
fig.savefig('figures/HYWIND_power_curve.pdf', bbox_inches='tight')

# Plot NREL power curve
WT = hes_off_object_oriented.process_models.WT(model='NREL')
fig, _ = WT.plot_power_curve()
fig.savefig('figures/NREL_power_curve.pdf', bbox_inches='tight')

# Plot yearly wind speed at Sleipner
filename = 'SLEIPNERWIND'
wind_data = hes_off_object_oriented.read_wind_data(filename)
fig, _ = hes_off_object_oriented.plot_wind_timeseries(wind_data, is_sorted=False)
fig.savefig('figures/sleipner_wind_timeseries.pdf', bbox_inches='tight')

# Plot yearly power generation at Sleipner
filename = 'SLEIPNERWIND'
wind_data = hes_off_object_oriented.read_wind_data(filename)
fig, _ = WT.plot_power_timeseries(wind_data["speed"], wind_data["time"], is_sorted=True)
fig.savefig('figures/sleipner_power_timeseries.pdf', bbox_inches='tight')

# Show the figure
plt.show()






