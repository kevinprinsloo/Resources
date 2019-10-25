"""
Parametric information and entropy estimation
"""

import numpy as np
import scipy as sp
import scipy.stats

def ctransform(x):
    "Rank and scale to [0, 1]"
    if len(x.shape) != 1:
        raise ValueError('Vector input required')

    xi = np.argsort(x)
    xr = np.argsort(xi)
    cx = (xr+1).astype(np.float) / (xr.size+1)

    return cx
    
    

def copnorm(x):
    """Copula-normalise a data set"""
    if len(x.shape) != 1:
        raise ValueError('Vector input required')

    cx = ctransform(x)
    cx = sp.stats.norm.ppf(cx)

    return cx

#def ent_g(x,biascorrect=True,demeaned=False):
    #"Entropy of a gaussian"

def info_gg(x,y,biascorrect=True):
    
    if len(x.shape)>1 or len(y.shape)>1:
        raise ValueError('Only vectors supported')
    
    Ntrl = len(x);
    if len(y) != Ntrl:
        raise ValueError('Different numbers of trials')

    C = np.cov(np.vstack((x,y)))

    Hx = 0.5*np.log(2*np.pi*np.e*C[0,0]) 
    Hy = 0.5*np.log(2*np.pi*np.e*C[1,1]) 
    Hxy = np.log(2*np.pi*np.e) + 0.5*np.log( C[0,0]*C[1,1] - C[0,1]**2 )
    I = Hx + Hy - Hxy
    return I


