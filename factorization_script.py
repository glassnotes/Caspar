#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# factorization_script.py: Factorize the matrix that has been inputted in the
#                          file user_matrix.py. 
#                                                                                  
# © 2017 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Caspar.
# Licensed under BSD-3-Clause                                                      
#                                                                                  

import sys
import numpy as np

from caspar import * 

# Read in input from file
from user_matrix import SUn_mat

np.set_printoptions(precision = 6)

# Get the dimension of the matrix.
n = SUn_mat.shape[0]

print("Original SU(" + str(n) + ") matrix is")
print(SUn_mat)

# Now perform the decomposition
parameters = sun_factorization(SUn_mat)

if parameters is not None:
    print("Factorization parameters: ")
    for param in parameters:
        print(param[0] + "\t" + str(param[1]))

# Reconstruct the matrix; set tiny numbers < 1e-15 to 0
# Uncomment if desired
"""
print("Reconstruct the matrix using output parameters: ")
recon = sun_reconstruction(n, parameters)
recon.real[np.abs(recon.real) < 1e-15] = 0.0
recon.imag[np.abs(recon.imag) < 1e-15] = 0.0
print(recon)
"""
