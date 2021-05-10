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
import numpy as np
import matplotlib.pyplot as plt
from hes_off_object_oriented import compute_specific_carbon_dioxide_emissions, convert_molar_to_mass_fraction

# Create output folder if it does not exist
if os.path.exists('figures/') is False:
    os.mkdir('figures')

# Define font settings
plt.rc('text', usetex=False)
plt.rcParams['font.family']      = 'serif'      # 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'
plt.rcParams['font.serif']       = 'cmr10'      # 'cmr10', 'palatino', 'times new roman'
plt.rcParams['mathtext.fontset'] = 'cm'         # 'cm' (latex style), 'stix' (times new roman style), 'stixsans'

# Define coal specifications
MW_coal  = np.asarray([12.00]) / 1e3
LHV_coal = np.asarray([32.80]) * 1e6
CR_coal  = np.asarray([1])
x_coal   = np.asarray([1.00])

# Define the diesel specifications (assuming that the average chemical formula is C12H23)
MW_diesel  = np.asarray([167.00]) / 1e3
LHV_diesel = np.asarray([45.6]) * 1e6
CR_diesel  = np.asarray([12])
x_diesel   = np.asarray([1.00])

# Define the natural gas specifications (composition order: CH4, C2H6, C3H8, C4H10)
MW_NG  = np.asarray([16.00, 30.00, 44.00, 59.00]) / 1e3
LHV_NG = np.asarray([55.50, 51.90, 50.35, 49.50]) * 1e6
CR_NG  = np.asarray([1, 2, 3, 4])
x_NG   = np.asarray([0.900, 0.050, 0.030, 0.020])
y_NG   = convert_molar_to_mass_fraction(x_NG, MW_NG)

# Compute the specific carbon dioxide emissions for different gas turbine efficiencies
efficiency = np.linspace(0.3, 1.00, 100)
specific_emissions_coal   = compute_specific_carbon_dioxide_emissions(efficiency, x_coal, MW_coal, LHV_coal, CR_coal)
specific_emissions_diesel = compute_specific_carbon_dioxide_emissions(efficiency, x_diesel, MW_diesel, LHV_diesel, CR_diesel)
specific_emissions_NG     = compute_specific_carbon_dioxide_emissions(efficiency, x_NG, MW_NG, LHV_NG, CR_NG)

# Plot the specific carbon dioxide emissions
fontsize = 14
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)
ax.set_xlabel('Conversion efficiency (%)', fontsize=fontsize, color='k', labelpad=fontsize)
ax.set_ylabel('Specific CO$_2$ emissions (g/kWh)', fontsize=fontsize, color='k', labelpad=fontsize)
# ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
# ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)

ax.plot(efficiency * 100, specific_emissions_coal, linewidth=1.25, linestyle='--', color='k',
        marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label='Coal')

ax.plot(efficiency * 100, specific_emissions_diesel, linewidth=1.25, linestyle='-.', color='k',
        marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label='Diesel')

ax.plot(efficiency * 100, specific_emissions_NG, linewidth=1.25, linestyle='-', color='k',
        marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label='Natural gas')

ax.set_ylim([-100, 1500])
ax.legend(ncol=1, loc='upper right', fontsize=fontsize-1, edgecolor='k', framealpha=1.0)
fig.savefig('figures/specific_carbon_dioxide_emissions_fossil.pdf', bbox_inches='tight')
fig.tight_layout()
