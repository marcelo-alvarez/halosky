import numpy as np
from numpy.random import *
from scipy.optimize import root_scalar
from scipy.interpolate import interp1d
from . import cosmology as co
import matplotlib.pyplot as plt

def fsky2fov_root(fov,fsky):
    return fsky - fov/4/np.pi*(np.cos((np.pi-fov)/2)-np.cos((np.pi+fov)/2))

def fsky2fov(fsky):
    return root_scalar(fsky2fov_root,args=[fsky],bracket=[0, np.pi], method='brentq').root

class lightcone:
    '''Lightcone'''
    def __init__(self, **kwargs):

        self.cosmo   = kwargs.get('cosmo',co.cosmology())
        self.zmin = kwargs.get('zmin',0.0)
        self.zmax = kwargs.get('zmax',5.0)
        self.Mmin = kwargs.get('Mmin',1e13)
        self.Mmax = kwargs.get('Mmax',1e16)
        self.fsky = kwargs.get('fsky',0.01)

        self.dz   = kwargs.get('dz', 0.1)

        if(self.fsky < 1.0):
            if(self.fsky > 0.5):
                raise Exception("field of view for fsky>0.5 not defined")
            self.fovrad = fsky2fov(self.fsky)
            self.fov    = 180. / np.pi * self.fovrad
            # center of fov is oriented along the x-axis (phi=0, theta=pi/2 in spherical polar coordinates)
            # negative phi corresponds to negative y, positive phi to positive y
            # theta < pi / 2 corresponds positive z, theta > pi / 2 to negative z
            self.phi0   = -self.fovrad / 2
            self.phi1   =  self.fovrad / 2
            self.theta0 = np.pi / 2 - self.fovrad / 2
            self.theta1 = np.pi / 2 + self.fovrad / 2

        else:
            self.phi0   = -np.pi
            self.phi1   =  np.pi
            self.theta0 = 0.0
            self.theta1 = np.pi

        self.clear()

    def clear(self):

        self.Nobj = 0
        self.M    = np.zeros(self.Nobj)
        self.x    = np.zeros(self.Nobj)
        self.y    = np.zeros(self.Nobj)
        self.z    = np.zeros(self.Nobj)

    def _getngtm(self,z,z0,z1,dndmfunc):

        co = self.cosmo

        Narr = 10000

        chi0  = co.chiofz(z0)
        chi1  = co.chiofz(z1)
        chi03 = chi0**3
        chi13 = chi1**3
        vol   = 4./3.*np.pi*(chi13-chi03)*self.fsky

        logMmin = np.log10(self.Mmin)-0.1
        logMmax = np.log10(self.Mmax)+0.1

        Marr    = np.logspace(logMmin,logMmax,Narr)
        dM      = (np.log(Marr[1])-np.log(Marr[0])) * Marr
        dnarr   = dndmfunc(Marr,z) * dM
        ngtmarr = np.cumsum(dnarr[::-1])[::-1] * vol # N(>M)

        ngtmfunc  = interp1d(Marr,ngtmarr)
        ngtmfunci = interp1d(ngtmarr,Marr)

        return ngtmfunc, ngtmfunci

    def populate(self, dndmfunc):

        print("\n populating catalog with halos")

        Mmin = self.Mmin
        Mmax = self.Mmax
        zmin = self.zmin
        zmax = self.zmax
        dz   = self.dz
        fsky = self.fsky
        co   = self.cosmo

        mu0  = np.cos(self.theta0)
        mu1  = np.cos(self.theta1)
        phi0 = self.phi0
        phi1 = self.phi1
        dmu  = mu1 - mu0
        dphi = phi1 - phi0

        Nz = int((zmax-zmin)/dz)

        # loop over redshift shells
        for iz in np.arange(Nz):

            z0 = zmin + iz * dz
            z1 = z0 + dz
            z  = z0 + dz / 2

            # cumulative mass function in shell
            ngtmfunc, ngtmfunci = self._getngtm(z,z0,z1,dndmfunc)

            # N(>Mmin in shell)
            Nbar = ngtmfunc(Mmin)

            # total number of halos in shell
            N = poisson(Nbar)

            if(N==0): return

            self.Nobj += N

            # mass of each halo sampled from distribution
            fN = uniform(size=N)
            M  = ngtmfunci(Nbar*fN)

            print("   z, Nshell, Ntot, Mmin, Mmax: ","{:4.2f}".format(z),"{:9d}".format(N),"{:9d}".format(self.Nobj),
                  "{:e}".format(M.min()),"{:e}".format(M.max()), end="\r", flush=True)

            # distance sampled in volume
            fV    = uniform(size=N)
            chi03 = co.chiofz(z0)**3
            chi13 = co.chiofz(z1)**3
            dchi3 = chi13 - chi03
            chi   = (chi03 + dchi3*fV)**(1./3.)

            # angles sampled over fov
            fmu   = uniform(size=N)
            fphi  = uniform(size=N)

            mu    = mu0 + fmu*dmu
            phi   = phi0 + fphi*dphi

            z = chi * mu
            r = chi * np.sqrt(1-mu**2)
            x = r * np.cos(phi)
            y = r * np.sin(phi)

            self.M = np.append(self.M, M).astype(np.float32)
            self.x = np.append(self.x, x).astype(np.float32)
            self.y = np.append(self.y, y).astype(np.float32)
            self.z = np.append(self.z, z).astype(np.float32)

        return

    def write_pksc(self,**kwargs):

        import pathlib
        import os

        catfile = kwargs.get('catfile','halos.pksc')

        M = self.M
        x = self.x
        y = self.y
        z = self.z

        co = self.cosmo

        Nobj = self.Nobj

        rho = 2.775e11 * co.omegam * co.h**2
        R = ((3.*M/4./np.pi/rho)**(1./3.)).astype(np.float32)

        header = np.asarray([self.Nobj, 0, 0]).astype(np.int32)

        path   = pathlib.Path(catfile).resolve()
        parent = path.parent
        parent.mkdir(parents=True, exist_ok=True)

        f = open(catfile,'wb')
        header.tofile(f)

        chunksize = 10000000
        end = 0
        remain = True
        print("\n writing",Nobj,"halos to file",os.path.basename(path))
        while remain:
            start = end
            end   = start + chunksize
            if end > Nobj:
                remain = False
                end = Nobj
            print("   start_index, end_index: ", start, end, end="\r", flush=True)
            outdata = np.column_stack((x[start:end],y[start:end],z[start:end],
                                       x[start:end],y[start:end],z[start:end],
                                       R[start:end],
                                       x[start:end],y[start:end],z[start:end]))
            outdata.tofile(f)
        print("\n",end="\r", flush=True)
