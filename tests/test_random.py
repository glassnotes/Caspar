#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# test_random.py: See how well the factorization works for Haar-random 
#                 unitaries. By "well", I mean that the parameters produce
#                 reconstructed matrices that have entries within some tolerance
#                 of the original matrix. 
#                                                                                  
# © 2017 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Caspar.                                         
# Licensed under BSD-3-Clause                                                      
#    

import numpy as np

from math import log2
from functools import reduce

from caspar import *

# Relevant parameters
TOL = 1e-8
NUM_TESTS = 100
MIN_N = 3 # Minimum and maximum dimension
MAX_N = 10 # Typically we start seeing some issues around 8, 9, 10 for dense mats

for n in range(MIN_N, MAX_N + 1):
    results = []

    for i in range(NUM_TESTS):
        # Generate a random SU(n) matrix.
        real_half = np.random.randn(n, n)
        comp_half = np.random.randn(n, n) 
        U = real_half + 1j * comp_half;  
        Q, R = np.linalg.qr(U)
        SUn_mat = Q
        SUn_mat[0,:] = SUn_mat[0,:] / np.linalg.det(SUn_mat) 

        params = sun_factorization(np.asmatrix(SUn_mat))
        recon = sun_reconstruction(n, params)

        if not np.allclose(SUn_mat, recon, atol = TOL):
            results.append(False)
        else:
            results.append(True)

    num_success = results.count(True) 
    
    print("Decomposition success for n = " + str(n) + " in " \
            + str(num_success) + "/" + str(NUM_TESTS) + " cases.")
