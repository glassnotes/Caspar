#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# matrixfile.py: Enter your matrix in this file as a numpy array.
#                To convert from the typical Mathematica format, it suffices
#                to replace the { } with [ ], complex variable I with *1j, and 
#                pi to np.pi.
#                                                                                  
# © 2017 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Caspar. 
# Licensed under BSD-3-Clause                                                      
#   

import numpy as np

# Modify me to whatever matrix you like!
# Currently n = 10 is about the largest Caspar can safely handle for very dense
# matrices, but it has been tested on up to 6-qubit Paulis (more sparse)
# successfully.

SUn_mat = np.matrix([[0, 0, 1.],
                    [np.exp(2 * 1j * np.pi / 3), 0., 0],
                    [0, np.exp(-2 * 1j * np.pi / 3), 0.]])
