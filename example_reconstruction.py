#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# example_reconstruction.py: Code snippet showing how an arbitrary set of 
#                            parameters and relevant modes can be used to 
#                            compute a unitary.
#                                                                                  
# © 2017 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Caspar.
# Licensed under BSD-3-Clause                                                      
#                                                                                  

import numpy as np

from sun_reconstruction import sun_reconstruction

np.set_printoptions(precision = 4)

# We can reconstruct arbitrary sequences as well by providing a list of tuples
# of parameters and modes on which they act. Modes are 1-indexed.
# Note that the order of multiplication is "as-written", i.e. (U_1, U_2, U_3) 
# will be multiplied together to form U_1 U_2 U_3, and *not* U_3 U_2 U_1.

# Here's some arbitrary parameters and modes they act on. 
n = 8
some_params = [("1,2", [1.23, 2.34, 0.999]), ("4,5", [0.3228328, 0.23324, -0.2228])]
random_transform = sun_reconstruction(n, some_params)

print("Reconstruction of matrix with some arbitrary parameters.")
print(random_transform)
