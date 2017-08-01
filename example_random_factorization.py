#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# example_random_factorization.py: A code snippets demonstrating our 
#                                  factorization scheme on a random matrix.
#                                                                                  
# © 2017 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Caspar.
# Licensed under BSD-3-Clause                                                      
#                                                                                  

import numpy as np

from sun_factorization import sun_factorization
from sun_reconstruction import sun_reconstruction

np.set_printoptions(precision = 4)

# Pick your dimension. Currently it works very well up until
# dimension 11, then there starts to be some issues with precision. 
n = 8

# Create a random SU(n) matrix; first find a random unitary
# and scale by its determinant to get something in SU(n).
real_half = np.random.randn(n, n)
comp_half = np.random.randn(n, n)
U = real_half + 1j * comp_half;
Q, R = np.linalg.qr(U)
SUn_mat = Q
SUn_mat[0,:] = SUn_mat[0,:]/np.linalg.det(SUn_mat)

print("Original random SU(" + str(n) + ") matrix is")
print(SUn_mat)

# Now perform the decomposition
parameters = sun_factorization(SUn_mat)

print("Factorization parameters: ")  
for param in parameters:
    print(param[0] + "\t" + str(param[1]))

#print("Matrix reconstructed using parameters")
#print(sun_reconstruction(n, parameters))
