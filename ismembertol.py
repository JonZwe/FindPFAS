# -*- coding: utf-8 -*-
"""
Function to find identicals within a tolerance
input: two arrays (vec_1 and vec_2) and absolute tolerance (tol)
output: booleans with positions of identical values 
(bool_1 for vec_1, bool_2 for vec_2)
vec_1[bool_1] = values in vec_1 that are identical with vec_2 within tol
"""
def ismembertolerance(vec_1,vec_2,tol):
    import numpy as np
    
    bool_1 = np.zeros(len(vec_1))
    for n in range(len(vec_1)):
        counter = 0
        for m in range(len(vec_2)):
            diff = abs(vec_1[n] - vec_2[m]) 
            if diff <= tol:
                counter = 1      
        bool_1[n] = counter
        
    bool_1 = bool_1 == 1
    
    # =================================
    
    bool_2 = np.zeros(len(vec_2))
    for n in range(len(vec_2)):
        counter = 0
        for m in range(len(vec_1)):
            diff = abs(vec_2[n] - vec_1[m]) 
            if diff <= tol:
                counter = 1      
        bool_2[n] = counter
        
    bool_2 = bool_2 == 1
    
    # =================================
    return  bool_1, bool_2   
    