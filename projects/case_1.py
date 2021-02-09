import numpy as np
import os
import sys

# BLADE_HOME = os.environ["BLADE_HOME"]
# sys.path.append('../HES-OFF/')
#
# from evaluate_sum import *

import hes_off as hf

a, b = 1, 2
c = hf.evaluate_sum(a,b)
print(c)

print(sys.path)

