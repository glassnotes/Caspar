#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# test_qubit_paulis.py: See how well the factorization works for randomly
#                       selected n-qubit Pauli matrices.
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
MIN_QUBITS = 2 # Test from min to max number of qubits
MAX_QUBITS = 6 # 6 starts to get pretty slow, even for 100 tests

I = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])
mapping = {0 : I, 1 : X, 2 : Y, 3 : Z}

for num_qubits in range(MIN_QUBITS, MAX_QUBITS + 1):
    results = []

    for i in range(NUM_TESTS):
        # Generate a random Pauli
        random_paulis = [mapping[x] for x in np.random.randint(0, 3, num_qubits)]

        U = np.asmatrix(reduce(np.kron, random_paulis))
        params = sun_factorization(U)
        recon = sun_reconstruction(2 ** num_qubits, params)

        if not np.allclose(U, recon, atol = TOL):
            results.append(False)
        else:
            results.append(True)

    num_success = results.count(True) 
    
    print("Decomposition success for " + str(num_qubits) + " qubits in " \
            + str(num_success) + "/" + str(NUM_TESTS) + " cases.")
