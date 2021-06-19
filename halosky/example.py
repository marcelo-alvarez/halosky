import halosky.hmf       as hshmf
import halosky.lightcone as hslc
import halosky.halopaint as hspaint
import halosky.cosmology as hscosmo

import numpy as np
import matplotlib.pyplot as plt

fsky=1.0
dz=0.1
zmin=0
zmax=4.5
mmin=1e13
mmax=1e16

print('hscosmo is ',hscosmo.h)
# create function dndm(M,z) with input mass in Msun, and return units 1/Mpc^3/Msun
# here we use Tinker et al. (2008) with Websky cosmology

#dndmofmz = hshmf.hmf.dndmofmz_tinker(mmin,mmax,zmin,zmax)

# create lightcone object
#lc = hslc.lightcone(fsky=fsky,zmin=zmin,zmax=zmax,dz=dz,Mmin=mmin,Mmax=mmax)

# sample halos from the halo mass function in the lightcone
#lc.populate(dndmofmz)

# write halos in pksc halo format used in Websky
#lc.write_pksc('rancat_test.pksc')
