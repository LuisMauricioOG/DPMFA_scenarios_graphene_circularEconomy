# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 21:50:43 2018

@author: Sana
"""

import numpy as np
import numpy.random as nr


#s=5
#
#[x,y,z] = [0.1, 0.2, 0.3]
#
#a = nr.triangular(x, y, y, s)
#b = nr.triangular(y, y, z, s)
#
##function = 'nr.triangular'
#
#c =np.concatenate([a, b])
#nr.choice(c)

#nr.seed(2250)

def TriangularBalanced(x,y,z,s1=None):
    s=10000
    Tright = nr.triangular(x, y, y, s)
    Tleft = nr.triangular(y, y, z, s)
    Tcomplete =np.concatenate([Tright, Tleft])
    return nr.choice(Tcomplete, s1)
    
#newFunction(0.1,0.2,0.3)
