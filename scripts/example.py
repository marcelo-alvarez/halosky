import halosky as hs
import numpy as np

def examinecat(catfile):
    f=open(catfile)
    N=np.fromfile(f,count=3,dtype=np.int32)[0]

    # only take first five entries for testing (there are ~8e8 halos total...)
    # comment the following line to read in all halos
    #N = 5

    print('reading halos...')
    catalog=np.fromfile(f,count=N*10,dtype=np.float32)
    print('reshaping ...')
    catalog=np.reshape(catalog,(N,10))

    x  = catalog[:,0];  y = catalog[:,1];  z  = catalog[:,2] # Mpc (comoving)
    R  = catalog[:,6] # Mpc

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
    hp.makemap(catfile=catfile,mapname=mapname,nproc=16,nodes=1)

    # clear halos from Lightcone
    lc.clear()
