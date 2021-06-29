import halosky as hs

hp = hs.halopaint.halopaint(nside=512)

catfile='/global/cfs/cdirs/mp107/exgal/websky/data/halos/websky_halos_10x10.pksc'

hp.makemap(catfile=catfile,mapname="websky-tsz",fov=10,npix=4096,nodes=1,nproc=16)

