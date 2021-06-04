import numpy as np
from numba import njit
from numba.core import types
from numba.typed import Dict

# Make array type.  Type-expression is not supported in jit
# functions.
float_array = types.float64[:]

@njit
def foo():
    # Make dictionary
    d = Dict.empty(
        key_type=types.unicode_type,
        value_type=float_array,
    )
    # Fill the dictionary
    d["posx"] = np.arange(3).astype(np.float64)
    d["posy"] = np.arange(3, 6).astype(np.float64)
    return d

d = foo()
# Print the dictionary
print(d)  # Out: {posx: [0. 1. 2.], posy: [3. 4. 5.]}