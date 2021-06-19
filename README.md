# Halosky
Creating halo based maps of the sky

## Installation
1. git clone https://github.com/marcelo-alvarez/halosky.git
2. cd halosky
3. pip install .

## Running
Example included here in [scripts/example.py](https://github.com/marcelo-alvarez/halosky/blob/master/scripts/example.py):
```
import halosky as hs
import numpy as np
import matplotlib.pyplot as plt

# create cosmology
cosmo = hs.cosmology.cosmology(
    omegab = 0.049,
    omegac = 0.261,
    h      = 0.68,
    ns     = 0.965,
    sigma8 = 0.81
)

# create hmf object
hmf = hs.hmf.hmf(cosmo=cosmo)

# create lightcone object
lc = hs.lightcone.lightcone(cosmo=cosmo,fsky=1.0,Mmin=5e14)

#create halopaint object
hp = hs.halopaint.halopaint(nside=1024)

Nran = 10
for i in range(Nran):
    # sample halos from the halo mass function in the lightcone
    lc.populate(hmf.dndmofmz)

    # write halos in pksc halo format used in Websky
    catfile = './catalogs/cat_'+f'{i:05d}'+'.pksc'
    lc.write_pksc(catfile=catfile)

    # make map
    mapname = './maps/tsz_'+f'{i:05d}'
    hp.makemap(catfile=catfile,mapname=mapname,n=1,N=1)

    # clear halos from Lightcone
    lc.clear()
```
