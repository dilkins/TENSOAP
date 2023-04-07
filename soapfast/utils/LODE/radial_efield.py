import sys,os
import numpy as np
import scipy.special as sc
import math
import time

def radial_efield(nmax,orthomatrix,sigma):
    """Compute external field contribution to local potential"""

    # compute radial integrals int_0^\infty dr r^3 R_n(r)
    radint = np.zeros(nmax)
    for n in range(nmax):
        inner = 0.5*sc.gamma(n+1.5)*(sigma[n]**2)**(n+1.5)
        radint[n] = 2**float(1.0+float(n)/2.0) * sigma[n]**(4+n) * sc.gamma(2.0+float(n)/2.0) / np.sqrt(inner)
   
    # orthogonalize radial integrals
    orthoradint = np.dot(orthomatrix,radint)

    return orthoradint
