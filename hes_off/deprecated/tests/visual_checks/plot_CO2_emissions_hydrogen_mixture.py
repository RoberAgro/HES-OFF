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

# Define the natural gas specifications (composition order: CH4, C2H6, C3H8, C4H10)
MW_NG  = np.asarray([16.00, 30.00, 44.00, 59.00]) / 1e3
LHV_NG = np.asarray([55.50, 51.90, 50.35, 49.50]) * 1e6
CR_NG  = np.asarray([1, 2, 3, 4])
x_NG   = np.asarray([0.900, 0.050, 0.030, 0.020])
y_NG   = convert_molar_to_mass_fraction(x_NG, MW_NG)

# Define hydrogen specifications
MW_H2  = 2/1e3
LHV_H2 = 120e6
CR_H2   = 0
molar_fraction_H2 = np.linspace(0, 1, 100)
mass_fraction_H2 = []
specific_emissions_H2 = []

# Define a range of conversion efficiencies
efficiency = np.asarray([0.30, 0.40, 0.50, 0.60, 0.70, 0.80])

for x_H2 in molar_fraction_H2:

    # Define the specifications of the fuel mixture
    LHV = np.append(LHV_NG, LHV_H2)
    MW  = np.append(MW_NG,  MW_H2)
    CR  = np.append(CR_NG,  CR_H2)
    x   = np.append((1-x_H2)*x_NG, x_H2)

    # Compute the mass fraction of hydrogen
    mass_fraction_H2.append(convert_molar_to_mass_fraction(x, MW)[-1])

    # Compute the specific carbon dioxide emissions of the fuel mixture
    specific_emissions_H2.append(compute_specific_carbon_dioxide_emissions(efficiency, x, MW, LHV, CR))


# Convert list of arrays into two-dimensional array
mass_fraction_H2 = np.asarray(mass_fraction_H2)
specific_emissions_H2 = np.asarray(specific_emissions_H2)


# Plot the specific carbon dioxide emissions as a function of the H2 molar fraction
fontsize = 14
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)
ax.set_xlabel('Molar fraction of hydrogen in the fuel (%)', fontsize=fontsize, color='k', labelpad=fontsize)
ax.set_ylabel('Specific CO$_2$ emissions (g/kWh)', fontsize=fontsize, color='k', labelpad=fontsize)
# ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
# ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
start_color  = np.asarray([0.10, 0.10, 0.10])
end_color    = np.asarray([0.85, 0.85, 0.85])

N = specific_emissions_H2.shape[-1]
for i in range(N):
    color = start_color + end_color * i / N
    label = ' ' + '{:.0f}'.format(efficiency[i] * 100) + '%'
    ax.plot(molar_fraction_H2 * 100, specific_emissions_H2[:, i], linewidth=1.25, linestyle='-', color=color,
            marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=label)

ax.set_ylim([0, 700])
ax.legend(ncol=2, loc='upper right', fontsize=fontsize-2, edgecolor='k', framealpha=1.0,handletextpad=0.1, columnspacing=1)
fig.tight_layout()
fig.savefig('figures/specific_carbon_dioxide_emissions_H2_molar.pdf', bbox_inches='tight')


# Plot the specific carbon dioxide emissions as a function of the H2 mass fraction
fontsize = 14
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)
ax.set_xlabel('Mass fraction of hydrogen in the fuel (%)', fontsize=fontsize, color='k', labelpad=fontsize)
ax.set_ylabel('Specific CO$_2$ emissions (g/kWh)', fontsize=fontsize, color='k', labelpad=fontsize)
# ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
# ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
for t in ax.xaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
for t in ax.yaxis.get_major_ticks(): t.label1.set_fontsize(fontsize)
start_color  = np.asarray([0.10, 0.10, 0.10])
end_color    = np.asarray([0.85, 0.85, 0.85])

N = specific_emissions_H2.shape[-1]
for i in range(N):
    color = start_color + end_color * i / N
    label = ' ' + '{:.0f}'.format(efficiency[i] * 100) + '%'
    ax.plot(mass_fraction_H2*100, specific_emissions_H2[:, i], linewidth=1.25, linestyle='-', color=color,
            marker=' ', markersize=4.5, markeredgewidth=1.25, markeredgecolor='k', markerfacecolor='w', label=label)

ax.set_ylim([0, 700])
ax.legend(ncol=2, loc='upper right', fontsize=fontsize-2, edgecolor='k', framealpha=1.0,handletextpad=0.1, columnspacing=1)
fig.tight_layout()
fig.savefig('figures/specific_carbon_dioxide_emissions_H2_mass.pdf', bbox_inches='tight')
plt.show()