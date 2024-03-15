
# This code tests the root_newton_raphson function

import unittest
import numpy as np
import sys
import math

from my_python_package.rootfind import root_newton_raphson

def f(x):
    return  (math.exp(-x)) - x

def f_prime(x):
    return  -1*(math.exp(-x)) - 1

x_initial = 0


class TestIntegrationNewton(unittest.TestCase):
     
    def test_newtonraphson(self):

        exact_root = 0.567143290
        exact_i = 4
        exact_error = [100/100 , 11.8/100, 0.147/100, 0.0000220/100, 0]

        computed_root, computed_i, computed_error  = root_newton_raphson(x_initial, f, f_prime)

        self.assertAlmostEqual(computed_root, exact_root, delta=1e1)
        self.assertAlmostEqual(computed_i, exact_i, delta=1e1)
        self.assertTrue(np.allclose(computed_error, exact_error, atol=1e-5))

       
if __name__ == '__main__':
    unittest.main()
