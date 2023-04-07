import sys
import numpy as np
import scipy.special as sc
import math
import time
from . import phasecomb,gcontra,fourier_integrals,fourier_integrals_MT2D,realspace_potential,potrefewald,potrefewaldmt2d

def fourier_ewald(sigewald,nat,nnmax,nspecies,lmax,centers,all_species,nneighmax,atom_indexes,cell,rcut,coords,all_radial,sigma,sg,nmax,orthomatrix,nside,iGx,imGx,Gval,Gvec,nG,MT2D):
    """return projections of the non-local field on basis functions"""

    volume = np.linalg.det(cell)

    start = time.time()

    start1 = time.time()
    # process coordinates 
    coordx = np.zeros((nat,nspecies,nat,3), dtype=float)
    nneigh = np.zeros((nat,nspecies),int)
    iat = 0
    ncentype = len(centers)
    # loop over species to center on
    for icentype in range(ncentype):
        centype = centers[icentype]
        # loop over centers of that species
        for icen in range(nneighmax[centype]):
            cen = atom_indexes[centype,icen]
            # loop over all the species to use as neighbours
            for ispe in range(nspecies):
                spe = all_species[ispe]
                # loop over neighbours of that species
                n = 0
                for ineigh in range(nneighmax[spe]):
                    neigh = atom_indexes[spe,ineigh]
                    coordx[iat,ispe,n,0] = coords[neigh,0] - coords[cen,0]
                    coordx[iat,ispe,n,1] = coords[neigh,1] - coords[cen,1]
                    coordx[iat,ispe,n,2] = coords[neigh,2] - coords[cen,2]
                    nneigh[iat,ispe] += 1
                    n += 1
            iat = iat + 1

    start2 = time.time()
    # combine phase factors
    phase = phasecomb.phasecomb(nat,nspecies,nneigh,nG,coordx,Gvec.T)
#    print "phases combinations :", time.time()-start2, "seconds"

    start3 = time.time()
    # compute analytic radial integrals and spherical harmonics
    alphaewald = 1.0/(2.0*sigewald**2)
    if MT2D:
        [orthoradint,harmonics] = fourier_integrals_MT2D.fourier_integrals_MT2D(nG,nmax,lmax,alphaewald,rcut,sigma,Gval,Gvec,orthomatrix) 
    else:
        [orthoradint,harmonics] = fourier_integrals.fourier_integrals(nG,nmax,lmax,alphaewald,rcut,sigma,Gval,Gvec,orthomatrix) 
    orthoradint = np.moveaxis(orthoradint, 0, -1)
#    print "fourier integrals :", time.time()-start3, "seconds"

    start4 = time.time()
    # perform contraction over G-vectors
    omega = gcontra.gcontra(nat,nspecies,nmax,lmax,nG,orthoradint,harmonics.T,2.0*phase)
#    print "G-contraction :", time.time()-start4, "seconds"

#    print "-----------------------------------------"
    print("reciprocal space potential computed in", time.time()-start, "seconds")

    omega *= 16.0*np.pi**2/volume

#    #TODO uncomment to recenter potential
#    if MT2D:
#        potatom = potrefewaldmt2d.potrefewaldmt2d(nG,Gval,Gvec.T,alphaewald)
#    else:
#        potatom = potrefewald.potrefewald(nG,Gval,Gvec.T,alphaewald)
#    potatom *= 4*np.pi/volume
#    potref = np.sqrt(2.0/np.pi)/sigewald
#    potdiff = nat * (potatom-potref)
# 
#    intradialbasis = np.zeros(nmax,float)
#    for n in range(nmax):
#        # normalization factor for primitive radial functions
#        normfact = np.sqrt(2.0/(sc.gamma(1.5+n)*sigma[n]**(3.0+2.0*n)))
#        intradialbasis[n] = normfact * 2.0**((1+n)/2.0) * sigma[n]**(3+n) * sc.gamma(0.5*(3.0+n))
#    omegaref = np.sqrt(4.0*np.pi) * np.dot(orthomatrix,intradialbasis) * potdiff
#    for n in range(nmax):
#        omega[:,:,n,0,0] -= omegaref[n]

#    # define 3D grid
#    nside = {}
#    spacing = 0.5
#    for nn in range(1,1000):
#        dx1 = cell[0,0]/nn
#        dx2 = cell[0,0]/(nn+1)
#        if dx2 <= spacing <= dx1:
#            nside[0] = nn
#        dy1 = cell[1,1]/nn
#        dy2 = cell[1,1]/(nn+1)
#        if dy2 <= spacing <= dy1:
#            nside[1] = nn
#        dz1 = cell[2,2]/nn
#        dz2 = cell[2,2]/(nn+1)
#        if dz2 <= spacing <= dz1:
#            nside[2] = nn
#    
#    npoints = 1
#    for i in range(3):
#        npoints *= nside[i]
#    dx = cell[0,0] / nside[0]  # bohr 
#    dy = cell[1,1] / nside[1]  # bohr 
#    dz = cell[2,2] / nside[2]  # bohr 
#
#    potential = realspace_potential.realspace_potential(nG,nat,nside,dx,dy,dz,alphaewald,Gval,Gvec,coords)
#
#    potential /= cell[0,0]*cell[1,1]*cell[2,2]
#
#    cubef = open("potential_screened.cube","w")
#    print("Reconstructed electron density",file=cubef)
#    print("CUBE FORMAT",file=cubef)
#    print(nat, 0.0, 0.0, 0.0,file=cubef)
#    metric = np.array([[-dx,0.0,0.0],[0.0,-dy,0.0],[0.0,0.0,-dz]])
#    for ix in range(3):
#        print(nside[ix], metric[ix,0], metric[ix,1], metric[ix,2],file=cubef)
#    for iat in range(nat):
#        print(79, float(79), coords[iat,0], coords[iat,1], coords[iat,2],file=cubef)
#    for igrid in range(npoints):
#        print(potential[igrid],file=cubef)
#    cubef.close()
#
#    print("all good so far!")
#    sys.exit(0)

    return omega
