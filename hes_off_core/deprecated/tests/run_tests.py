#!/usr/bin/env python3
import pytest

# Define the list of tests
tests_process   = ["test_process_models.py"]
tests_grid      = ["test_grid_models.py"]
test_combustion = ["test_models.py"]
test_other      = ["test_utilities.py"]
tests_list  = tests_process + tests_grid + test_combustion + test_other

# Run pytest when this script is executed
pytest.main(tests_list)
