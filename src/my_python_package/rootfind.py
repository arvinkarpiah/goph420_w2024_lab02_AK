# This function computes the root of a function g using the Newton-Raphson method
# Inputs : initial guess, function and derivative of function
# Outputs: roots, number of iterations to converge and approximate relative error for each iteration

import numpy as np
import math

def root_newton_raphson(x0, f, dfdx):
    x = np.zeros(10000, dtype=float)
    x[0] = x0
    epsilon_x = np.zeros(10000, dtype=float)
    epsilon_x[0] = 1
    epsilon_s = 0.5 * 10**(-4)

    i = 0  
    while epsilon_x[i] > epsilon_s:
            
                x[i+1] = x[i] - (f(x[i]) / dfdx(x[i]))
                epsilon_x[i+1] = abs((x[i+1] - x[i]) / (x[i+1]))
                i += 1
     
    return float(x[i]), int(i), epsilon_x[1:i+2]  
