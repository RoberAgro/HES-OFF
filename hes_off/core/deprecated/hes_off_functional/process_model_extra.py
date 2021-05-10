import pdb
from .process_model import *


# @nb.jit(nopython=True, cache=True)
def evaluate_process_model(HEAT_DEMAND, POWER_DEMAND,
                           GT_MODEL, GT_UNITS, GT_MAX_H2,
                           WT_MODEL, WT_RATED_POWER, WT_REF_HEIGHT, WT_HUB_HEIGHT,
                           EL_MODEL, EL_RATED_POWER, EL_EFFICIENCY,
                           FC_MODEL, FC_RATED_POWER, FC_EFFICIENCY,
                           H2_CAPACITY, H2_INITIAL_FRACTION, H2_RECHARGE_THRESHOLD, H2_COFIRE_THRESHOLD,
                           WIND_SPEED, WIND_TIME,
                           natural_gas, hydrogen):

    #  Check the size of the power demand and heat demand arrays
    POWER_DEMAND, HEAT_DEMAND = np.atleast_1d(POWER_DEMAND), np.atleast_1d(HEAT_DEMAND)
    if POWER_DEMAND.size != HEAT_DEMAND.size:
        raise Exception("The number of elements of POWER_DEMAND and HEAT_DEMAND must be the same")

    # Initialize arrays to store the solution at each instance of the year
    p, n = POWER_DEMAND.size, WIND_SPEED.size
    flag = np.empty((p,n))                   # Type of operational strategy at each time instance
    GT_power = np.empty((p,n))               # Power generated in the gas turbines (W)
    WT_power = np.empty((p,n))               # Power generated in the wind turbines (W)
    FC_power = np.empty((p,n))               # Power generated in the fuel cell system (W)
    EL_power = np.empty((p,n))               # Power consumed in the electrolyzer system (W)
    H2_level = np.empty((p,n))               # Mass of hydrogen in the storage system (kg)
    H2_produced = np.empty((p,n))            # Mass of hydrogen produced in the electrolyzer (kg)
    H2_utilized = np.empty((p,n))            # Mass of hydrogen utilized in the fuel cell + gas turbine (kg)
    CO2_produced = np.empty((p,n))           # Mass of carbon dioxide emitted to the atmosphere (kg)
    power_deficit = np.empty((p,n))          # Power demand not satisfied (W)
    # energy_surplus = np.empty((p,n))         # Extra wind energy that is dissipated (J)

    # Initialize time array (must have the same shape as the other arrays to return within dictionary)
    times = np.empty((p,n),dtype=np.int32)
    for i in range(p):
        times[i,:] = np.arange(0, n)

    # Create the natural gas + hydrogen mixture used in the gas turbines
    blend_NG_H2 = create_fluid_mixture(fluids=(natural_gas, hydrogen),
                                       fractions=(1.00 - GT_MAX_H2, GT_MAX_H2),
                                       fraction_type="molar")

    # Loop over the time periods
    for p, (power_demand, heat_demand) in enumerate(zip(POWER_DEMAND, HEAT_DEMAND)):

        # Compute the initial level of hydrogen in the storage system
        H2_level[p, 0] = H2_CAPACITY * H2_INITIAL_FRACTION

        # Compute the minimum GT load required to satisfy the heat demand
        GT_power_min = compute_GT_power_from_heat(model=GT_MODEL, number_of_units=GT_UNITS, heat_output=heat_demand)[0]

        # Compute the maximum GT load of the current GT model
        GT_power_max = compute_GT_maximum_power(model=GT_MODEL, number_of_units=GT_UNITS)

        # Compute wind power over the year
        WT_power_available = compute_WT_power_output(model=WT_MODEL, hub_height=WT_HUB_HEIGHT, ref_height=WT_REF_HEIGHT, rated_power=WT_RATED_POWER, wind_speed=WIND_SPEED)

        for t in times[p]:

            # Use a run-out-of-steam strategy to supply the power demand
            if GT_power_min + WT_power_available[t] >= power_demand:

                # Case 1: Use GT_min and WT to satisfy the power demand (use EL to recharge H2)
                if H2_level[p,t] < H2_CAPACITY:
                    flag_current = 1
                    GT_power_current = GT_power_min
                    EL_power_current = np.minimum(EL_RATED_POWER, GT_power_min + WT_power_available[t] - power_demand)
                    WT_power_current = power_demand + EL_power_current - GT_power_current
                    FC_power_current = 0.00

                # Case 2: Use GT_min and WT to satisfy the power demand (do not use EL to recharge H2)
                else:
                    flag_current = 2
                    GT_power_current = GT_power_min
                    WT_power_current = power_demand - GT_power_current
                    EL_power_current = 0.00
                    FC_power_current = 0.00

            elif GT_power_max + WT_power_available[t] >= power_demand:

                # Case 3: Use GT and WT to satisfy the power demand (use GT+EL to recharge H2)
                if H2_level[p,t] < H2_RECHARGE_THRESHOLD * H2_CAPACITY:
                    flag_current = 3
                    WT_power_current = WT_power_available[t]
                    EL_power_current = np.minimum(EL_RATED_POWER, GT_power_max + WT_power_current - power_demand)
                    GT_power_current = np.minimum(GT_power_max, power_demand + EL_power_current - WT_power_current)
                    FC_power_current = 0.00

                # Case 4: Use GT and WT to satisfy the power demand (do not use GT+EL to recharge H2)
                else:
                    flag_current = 4
                    WT_power_current = WT_power_available[t]
                    GT_power_current = power_demand - WT_power_current
                    EL_power_current = 0.00
                    FC_power_current = 0.00

            else:

                # Case 5: Use GT, WT and FC to satisfy the power demand (there is hydrogen available)
                if H2_level[p,t] > 0.00:
                    flag_current = 5
                    WT_power_current = WT_power_available[t]
                    GT_power_current = GT_power_max
                    FC_power_current = power_demand - WT_power_current - GT_power_current
                    EL_power_current = 0.00

                # Case 6: The GT and WT cannot satisfy the power demand (there is no hydrogen available)
                else:
                    flag_current = 6
                    WT_power_current = WT_power_available[t]
                    GT_power_current = GT_power_max
                    FC_power_current = 0.00
                    EL_power_current = 0.00


            # Determine the type of fuel used in the gas turbines
            use_fuel_blend = H2_level[p,t] > H2_COFIRE_THRESHOLD * H2_CAPACITY
            if use_fuel_blend:
                GT_fuel = blend_NG_H2
                H2_mass_fraction = GT_fuel["y"][-1]     # Hydrogen is the last component
            else:
                GT_fuel = natural_gas
                H2_mass_fraction = 0.00

            # Compute gas turbine efficiency
            GT_efficiency = compute_GT_efficiency(model=GT_MODEL, number_of_units=GT_UNITS, power_output=GT_power_current)

            # Compute the carbon dioxide emissions (kg/s)
            CO2_produced_current = compute_GT_carbon_dioxide_emissions(power_output=GT_power_current, conversion_efficiency=GT_efficiency, fuel=GT_fuel)[0]*3600

            # Compute the mass flow rate of hydrogen co-fired in the gas turbine (kg/s)
            H2_cofired_current = compute_GT_fuel_consumption(power_output=GT_power_current, conversion_efficiency=GT_efficiency, fuel=GT_fuel)[0]*H2_mass_fraction

            # Compute mass flow rate of hydrogen fed to the fuel cell system (kg/s)
            H2_utilized_current = compute_FC_hydrogen_consumption(model=FC_MODEL, efficiency_coefficients=FC_EFFICIENCY, rated_power=FC_RATED_POWER, power_output=FC_power_current)[0]

            # Compute the mass flow rate of hydrogen produced in the electrolyzer system (kg/s)
            H2_produced_current = compute_EL_hydrogen_production(model=EL_MODEL, efficiency_coefficients=EL_EFFICIENCY, rated_power=EL_RATED_POWER, power_input=EL_power_current)[0]

            # Evaluate the power balance (W)
            # power_deficit_current = np.maximum(0.0, power_demand + EL_power_current - WT_power_current - GT_power_current - FC_power_current)
            power_deficit_current = power_demand + EL_power_current - WT_power_current - GT_power_current - FC_power_current

            # Compute the hydrogen level for the next time instance (skip last time step computation)
            if t < times[p,-1]:
                H2_level[p,t+1] = H2_level[p,t] + (H2_produced_current-H2_utilized_current-H2_cofired_current) * 3600

            # Store the current solution in its corresponding array
            flag[p,t] = flag_current
            GT_power[p,t] = GT_power_current
            WT_power[p,t] = WT_power_current
            FC_power[p,t] = FC_power_current
            EL_power[p,t] = EL_power_current
            power_deficit[p,t] = power_deficit_current
            CO2_produced[p,t] = CO2_produced_current
            H2_produced[p,t] = H2_produced_current
            H2_utilized[p,t] = H2_utilized_current + H2_cofired_current

    # Store the results in a dictionary
    result_dict = {"flag":           flag*1.00,     # Conversion from integer to float for Numba
                   "times":          times*1.00,    # Conversion from integer to float for Numba
                   "GT_power":       GT_power,
                   "WT_power":       WT_power,
                   "FC_power":       FC_power,
                   "EL_power":       EL_power,
                   "H2_produced":    H2_produced,
                   "H2_utilized":    H2_utilized,
                   "H2_level":       H2_level,
                   "CO2_produced":   CO2_produced,
                   "power_deficit":  power_deficit}

    return result_dict


