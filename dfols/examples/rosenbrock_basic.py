# DFO-LS example: minimize the Rosenbrock function
from __future__ import print_function
import numpy as np
import dfols

# Define the objective function
def rosenbrock(x):
    return np.array([10.0 * (x[1] - x[0] ** 2), 1.0 - x[0]])

# Define the starting point
x0 = np.array([-1.2, 1.0])

# Set random seed (for reproducibility)
np.random.seed(0)

# For optional extra output details
import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

# User params
userParams = {'logging.save_poisedness': True, 'logging.save_diagnostic_info': True, 'logging.save_xk': True, 'logging.save_rk': True, 'logging.n_to_print_whole_x_vector': 3}

# Call DFO-LS
soln = dfols.solve(rosenbrock, x0, user_params=userParams)
#soln = dfols.solve(rosenbrock, x0)

# Display output
print(soln)
print(soln.diagnostic_info)
soln.diagnostic_info.to_csv('diagnostics.txt')
