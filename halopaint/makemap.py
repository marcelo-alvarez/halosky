#!/usr/bin/env python

import subprocess
import sys

path = './'

pkscfile = path + 'catalogs/websky_halos.pksc'
tabfile  = path + 'tables/proftab_deltac_planck.bin'
exe      = path + 'halopaint-f/bin/pks2map'

n=64
N=4
nside=4096
npix=1024
npixc=4096
fov=3 # field of view of map in degrees
zmin=0.0
zmax=5
mmin=1e12
scramble=0
virflag=2

#for field in ['ksz','tsz','kap']:
for field in ['tsz']:
    mappath=path+'maps/'+field+'_p2m'
    runargs=['srun','-n',str(n),'-N',str(N),'--account','mp107',exe,
         pkscfile,mappath,tabfile,field,str(zmin),str(zmax),str(mmin),str(nside),
             '1',str(scramble),str(npix),str(fov),'0',str(npixc),str(virflag)]
    fitsfile=mappath+'_hp.fits'
    print(runargs)
    subprocess.run(['rm','-f',fitsfile]) # cfitsio does not like overwriting fits files
    subprocess.run(runargs)
