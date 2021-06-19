import subprocess
import sys
import pathlib

abpath     = lambda filename: pathlib.Path(filename).resolve()
parentpath = lambda filename: pathlib.Path(filename).resolve().parent

def mkdir(filename):
    path = parentpath(filename)
    path.mkdir(parents=True, exist_ok=True)

class halopaint:
    '''Halo painting'''

    def __init__(self, **kwargs):

        self.nside = kwargs.get('nside', 4096)

        self.tablepath = kwargs.get('tablepath','/global/cfs/cdirs/mp107/exgal/software/halopaint-f/tables')
        self.binpath   = kwargs.get('binbath',  '/global/cfs/cdirs/mp107/exgal/software/halopaint-f/bin')

        self.exe     = self.binpath   + '/' + 'pks2map'
        self.tabfile = self.tablepath + '/' + 'proftab_deltac_planck.bin'

    def makemap(self,**kwargs):

        import os
        
        catfile = kwargs.get('catfile','halos.pksc')
        mapname = kwargs.get('mapname','tsz')

        mkdir(mapname)

        n        = kwargs.get('n',          64)
        N        = kwargs.get('N',           4)
        npix     = kwargs.get('npix',     1024)
        npixc    = kwargs.get('npixc',    4096)
        fov      = kwargs.get('fov',         3) # field of view of map in degrees
        zmin     = kwargs.get('zmin',      0.0)
        zmax     = kwargs.get('zmax',        5)
        mmin     = kwargs.get('mmin',     1e12)
        scramble = kwargs.get('scramble',    0)
        virflag  = kwargs.get('virflag',     2)

        exe     = self.exe
        tabfile = self.tabfile

        nside = self.nside

        tmpname='tsz'
        runargs=['srun','-n',str(n),'-N',str(N),exe,
                 catfile,tmpname,tabfile,'tsz',str(zmin),str(zmax),str(mmin),str(nside),
                     '1',str(scramble),str(npix),str(fov),'0',str(npixc),str(virflag)]
        fitsfile=tmpname+'_hp.fits'
        outfile=mapname+'.fits'
        logname=mapname+'.log'

        logfile = open(logname,'w')

        print(' making map from catalog...')
        subprocess.run(['rm','-f',fitsfile]) # cfitsio does not like overwriting fits files
        subprocess.run(runargs,stdout=logfile)
        subprocess.run(['mv',fitsfile,outfile])
        print(' map written to ',os.path.basename(outfile))
