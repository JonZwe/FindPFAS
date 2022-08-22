"""
Function to find identicals within a tolerance
input: two arrays (vec_1 and vec_2) and absolute tolerance (tol)
output: booleans with indices of identical values 
e.g.: vec_1[bool_1] = values in vec_1 that are identical with vec_2 within tol
"""
def ismembertolerance(vec_1,vec_2,tol):
    import numpy as np
    
    bool_1 = np.zeros((len(vec_2), len(vec_1)))
    for n in range(len(vec_2)):
        bool_1[n] =  abs(vec_2[n] - vec_1) <= tol
    bool_1 = np.sum(bool_1, axis = 0).astype('bool')
    
    
    bool_2 = np.zeros((len(vec_1), len(vec_2)))
    for n in range(len(vec_1)):
        bool_2[n] =  abs(vec_1[n] - vec_2) <= tol
    bool_2 = np.sum(bool_2, axis = 0).astype('bool')


    return  bool_1, bool_2
