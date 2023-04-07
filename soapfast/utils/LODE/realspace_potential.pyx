cimport cython # the main Cython stuff

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cython.parallel cimport prange

cimport scipy.special.cython_special as csc # Cython interfaces to Scipy Special functions (sph_harm, gamma, hyp1f1)

from libc.math cimport sqrt,abs
from libc.stdio cimport printf

cdef extern void zero_array (double *array, double value, size_t n)

# NumPy and Scipy Special Functions
import sys
import numpy as np
import scipy.special as sc

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _realspacepotential(long nG,
                              long nat,
                              long nx,
                              long ny,
                              long nz,
                              double dx,
                              double dy,
                              double dz,
                              double alpha,
                              double[:] Gval,
                              double[:,:] Gvec,
                              double[:,:] coords,
                              double[:] potential 
                ):

    # Py_ssize_t is the correct C type for Python array indexes
    cdef Py_ssize_t iG, iat, ix, iy, iz, igrid
    cdef double G2, Lz, Gs, x, y, z, arg

    Lz = abs(2*np.pi/Gvec[1,2])
  
    for iG in xrange(1,nG):
        G2 = Gval[iG]**2
        Gs = np.sqrt(Gvec[iG,0]**2+Gvec[iG,1]**2)
        fourierpot = np.exp(-G2/(4.0*alpha))/G2
        #screeningpotz = np.exp(-G2/(4.0*alpha))*np.real(sc.erfc((alpha*Lz+1.0j*Gvec[iG,2])/(2*np.sqrt(alpha)))) + np.cos(Gvec[iG,2]*Lz/2.0) * ( np.exp(-Gs*Lz/2.0) - 0.5*np.exp(-Gs*Lz/2.0)*sc.erfc((alpha*Lz-Gs)/(2*np.sqrt(alpha))) - 0.5*np.exp(+Gs*Lz/2.0)*sc.erfc((alpha*Lz+Gs)/(2*np.sqrt(alpha))) )
        #fourierpot -= screeningpotz/G2
        igrid = 0
        for ix in xrange(nx):
            x = ix*dx
            for iy in xrange(ny):
                y = iy*dy
                for iz in xrange(nz):
                    z = iz*dz
                    for iat in xrange(1):
                        arg = Gvec[iG,0]*(x-coords[iat,0])+Gvec[iG,1]*(y-coords[iat,1])+Gvec[iG,2]*(z-coords[iat,2])
                        potential[igrid] = potential[igrid] + 2.0 * fourierpot * np.cos(arg)
                    igrid = igrid + 1  

#---------------------------------------------------------------------------------------------

def realspace_potential(nG,nat,nside,dx,dy,dz,alpha,Gval,Gvec,coords):
    """Compute long-range part of the potential in real space."""
    nx = nside[0]
    ny = nside[1]
    nz = nside[2]
    potential = np.zeros(nx*ny*nz,float)
    _realspacepotential(nG,nat,nx,ny,nz,dx,dy,dz,alpha,Gval,Gvec,coords,potential)
    return potential 
