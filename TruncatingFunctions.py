# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 10:03:30 2018

@author: dew
"""
import numpy as np
    
### FUNCTION FOR TRUNCATING TRAPEZOIDAL DISTRIBUTIONS
def TrapezTrunc(TC1, TC2, spread1, spread2, N, linf, lsup):
    from scipy.stats import trapz 
        
    # define variables for trapezoidal distribution
    A = TC1*(1-spread1)
    B = TC2*(1+spread2)
    
    c = (TC1-A)/(B-A)
    d = (TC2-A)/(B-A)
    loc = A
    scale = B-A
    
    dist = trapz.rvs(c, d, loc, scale, N)
    
    truncdist = [i for i in dist if i >= linf]
    truncdist = [i for i in truncdist if i <= lsup]
    
    
    while len(truncdist) < N:
         adddist = trapz.rvs(c, d, loc, scale, N-len(truncdist))
         truncadddist = [i for i in adddist if i >= linf]
         truncadddist = [i for i in truncadddist if i <= lsup]
         
         truncdist = truncdist + truncadddist
    
    return np.asarray(truncdist)


### FUNCTION FOR TRUNCATING TRIANGULAR DISTRIBUTIONS
def TriangTrunc(TC1, spread1, N, linf, lsup):
    import numpy.random as nr
            
    # define variables for trapezoidal distribution
    A = TC1*(1-spread1)
    B = TC1*(1+spread1)
    
    dist = nr.triangular(A, TC1, B, N)
    
    # remove all that's not in the proper range - we end with a distribution with a length < N
    truncdist = [i for i in dist if i >= linf]
    truncdist = [i for i in truncdist if i <= lsup]
    
    
    # Create the "same" triangular distribution with the missing number of samples, and remove all that's not in the range.
    # "while" -> do that until you have N samples. 
    while len(truncdist) < N:
         adddist = nr.triangular(A, TC1, B, N-len(truncdist))
         truncadddist = [i for i in adddist if i >= linf]
         truncadddist = [i for i in truncadddist if i <= lsup]
         
         # concatenate both triangular distributions
         truncdist = truncdist + truncadddist
    
    return np.asarray(truncdist)
	
    
def TriangTruncDet(min, TC1, max, N, linf, lsup):
	import numpy.random as nr
	import numpy as np
	
	dist = nr.triangular(min, TC1, max, N)
	
	truncdist = [i for i in dist if i >= linf]
	truncdist = [i for i in truncdist if i <= lsup]
	
	while len(truncdist) < N:
		adddist = nr.triangular(min, TC1, max, N-len(truncdist))
		truncadddist = [i for i in adddist if i >= linf]
		truncadddist = [i for i in truncadddist if i <= lsup]
		
		truncdist = truncdist + truncadddist
		
	return np.asarray(truncdist)

# for testing
#import matplotlib.pyplot as plt
#plt.hist(TrapezTrunc(TC1=0.2, TC2=0.5, spread1=1.5, spread2=0.4, N=100, linf=0, lsup=1), bins=50)