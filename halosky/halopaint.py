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

        self.tablepath = kwargs.get('tablepath','/global/cfs/cdirs/mp107/exgal/users/malvarez/halopaint-f/tables')
        self.binpath   = kwargs.get('binbath',  '/global/cfs/cdirs/mp107/exgal/users/malvarez/halopaint-f/bin')

        self.exe     = self.binpath   + '/' + 'pks2map'
        self.tabfile = self.tablepath + '/' + 'proftab_websky.bin'

    def loadcutout(self,filename):
        import array
        import math
        import numpy
        import sys

        mapfile = open(filename,"rb")

        n=array.array('i')
        fov=array.array('f')
        vals=array.array('f')
        
        n.fromfile(mapfile,2)
        fov.fromfile(mapfile,2)
        
        print('dimensions of map are:  ',n[0],'x',n[1])
        print('field of view is:       ',fov[0]/2./math.pi*360.,'x',
              fov[1]/2./math.pi*360.,' degrees')

        vals.fromfile(mapfile,n[0]*n[1])

        vals=numpy.reshape(vals,(n[0],n[1]))

        return vals

    def makemap(self,**kwargs):

        import os

        catfile = kwargs.get('catfile','halos.pksc')
        mapname = kwargs.get('mapname','tsz')

        mkdir(mapname)

        n        = kwargs.get('nproc',      64)
        N        = kwargs.get('nodes',       4)
        npix     = kwargs.get('npix',     1024)
        npixc    = kwargs.get('npixc',    4096)
        fov      = kwargs.get('fov',         3) # field of view of map in degrees
        zmin     = kwargs.get('zmin',      0.0)
        zmax     = kwargs.get('zmax',        5)
        mmin     = kwargs.get('mmin',     1e12)
        scramble = kwargs.get('scramble',    0)
        virflag  = kwargs.get('virflag',     2)
        parallel = kwargs.get('parallel', True)

        exe     = self.exe
        tabfile = self.tabfile
        nside = self.nside

        tmpname='tsz'
        runargs=[exe,catfile,tmpname,tabfile,'tsz',str(zmin),str(zmax),str(mmin),
                 str(nside),'1',str(scramble),str(npix),str(fov),'0',str(npixc),
                 str(virflag)]
        if parallel:
            runargs = ['srun','-n',str(n),'-N',str(N)] + runargs
        print(runargs)
        fitsfile=tmpname+'_hp.fits'
        mapfile=tmpname+'_fs.map'
        dpfile=tmpname+'_dp.bin'
        carfile=tmpname+'_car.map'
        
        outfile_hp=mapname+'.fits'
        outfile_fs=mapname+'.map'
        
        logname=mapname+'.log'

        logfile = open(logname,'w')

        print(' making map from catalog...')
        subprocess.run(['rm','-f',fitsfile]) # cfitsio does not like overwriting fits files
        subprocess.run(runargs,stdout=logfile)
        subprocess.run(['mv',fitsfile,outfile_hp])
        subprocess.run(['mv',mapfile,outfile_fs])
        subprocess.run(['rm','-f',dpfile,carfile]) # remove unused files
        print(' maps written to ',os.path.basename(outfile_hp),' and ',
              os.path.basename(outfile_fs))
