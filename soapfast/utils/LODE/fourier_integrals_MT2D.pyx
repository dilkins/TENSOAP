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
cdef void _fintegrals(long nG,
                  long nmax,
                  long lmax,
                  double alpha,
                  double rcut,
                  double[:] sigma,
                  double[:] Gval,
                  double[:,:] Gvec,
                  double[:,:] radint,
                  double[:] intradialbasis,
                  double[:,:] orthomatrix, 
                  double[:,:] prefacts, 
                  double[:,:,:] orthoradint,
                  complex[:,:] harmonics 
                ):

    # Py_ssize_t is the correct C type for Python array indexes
    cdef Py_ssize_t iG, n, l, lm, im, n1, n2
    cdef double fourierpot, G2, normfact, th, ph , Lz, Gs, besselaverage, besselarg

    for n in xrange(nmax):
        # normalization factor for primitive radial functions
        normfact = np.sqrt(2.0/(sc.gamma(1.5+n)*sigma[n]**(3.0+2.0*n)))
        intradialbasis[n] = normfact * 2.0**((1+n)/2.0) * sigma[n]**(3+n) * sc.gamma(0.5*(3.0+n)) 
        for l in xrange(lmax+1):
            # precompute common prefactors for each (n,l) pair
            prefacts[n,l] =  normfact * np.sqrt(np.pi) * 2.0**((n-l-1)/2.0) * sigma[n]**(3.0+l+n) * sc.gamma(0.5*(3.0+l+n)) / sc.gamma(1.5+l)

    Lz = abs(2*np.pi/Gvec[1,2])
  
    for iG in xrange(1,nG):
        G2 = Gval[iG]**2
        Gs = np.sqrt(Gvec[iG,0]**2+Gvec[iG,1]**2)
        fourierpot = np.exp(-G2/(4.0*alpha))/G2
        screeningpotz = np.exp(-G2/(4.0*alpha))*np.real(sc.erfc((alpha*Lz+1.0j*Gvec[iG,2])/(2*np.sqrt(alpha)))) + np.cos(Gvec[iG,2]*Lz/2.0) * ( np.exp(-Gs*Lz/2.0) - 0.5*np.exp(-Gs*Lz/2.0)*sc.erfc((alpha*Lz-Gs)/(2*np.sqrt(alpha))) - 0.5*np.exp(+Gs*Lz/2.0)*sc.erfc((alpha*Lz+Gs)/(2*np.sqrt(alpha))) )
        fourierpot -= screeningpotz/G2
        for n in xrange(nmax):
            arghyper = -0.5*G2*(sigma[n]**2)
            for l in xrange(lmax+1):
                # radial integral
                radint[l,n] = prefacts[n,l] * sc.hyp1f1(0.5*(3.0+l+n), 1.5+l, arghyper)
                radint[l,n] *= fourierpot * Gval[iG]**l 
        # compute polar angles at G-vectors directions
        th = np.arccos(Gvec[iG,2]/Gval[iG])
        ph = np.arctan2(Gvec[iG,1],Gvec[iG,0])
        lm = 0
        for l in xrange(lmax+1):
            for n1 in xrange(nmax):
                for n2 in xrange(nmax):
                    # orthogonalize radial integrals with Loewdin
                    orthoradint[iG,l,n1] = orthoradint[iG,l,n1] + orthomatrix[n1,n2]*radint[l,n2]
            for im in xrange(2*l+1):
                # compute spherical harmonics at G-vector directions
                harmonics[iG,lm] = np.conj(sc.sph_harm(im-l,l,ph,th)) * 1.0j**l
                lm += 1


#-------------------------------------------------------------------------------------------------------------------------------------------

def fourier_integrals_MT2D(nG,nmax,lmax,alpha,rcut,sigma,Gval,Gvec,orthomatrix):
    """return radial integrals"""
    radint = np.zeros((lmax+1,nmax),float)
    orthoradint = np.zeros((nG,lmax+1,nmax),float)
    prefacts = np.zeros((nmax,lmax+1),float)
    intradialbasis = np.zeros(nmax,float)
    harmonics = np.zeros((nG,(lmax+1)**2), dtype=complex)
    _fintegrals(nG,nmax,lmax,alpha,rcut,sigma,Gval,Gvec,radint,intradialbasis,orthomatrix,prefacts,orthoradint,harmonics)
    return [orthoradint,harmonics]
