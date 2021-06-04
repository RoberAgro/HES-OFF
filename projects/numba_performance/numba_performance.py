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
import hes_off_core
import hes_off_core.deprecated.hes_off_functional as hes_off_functional
import hes_off_core.deprecated.hes_off_object_oriented as hes_off_object_oriented

# Evaluate the process model (object oriented)
IN = hes_off_core.get_defaults()
EnergySystem = hes_off_object_oriented.IntegratedModel(IN)
for i in range(5):
    with hes_off_core.Timer() as timer:
        EnergySystem.evaluate_process_model()

# Evaluate the process model (functional without Numba)
IN = hes_off_core.get_defaults()
EnergySystem = hes_off_functional.IntegratedModel(IN)
for i in range(5):
    with hes_off_core.Timer() as timer:
        EnergySystem.evaluate_process_model()

# Evaluate the process model (functional with Numba)
IN = hes_off_core.get_defaults()
EnergySystem = hes_off_core.IntegratedModel(IN)
for i in range(5):
    with hes_off_core.Timer() as timer:
        EnergySystem.evaluate_process_model()
