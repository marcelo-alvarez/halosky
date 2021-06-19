from scipy.interpolate import *
import numpy as np

class cosmology:
    '''Cosmology'''

    def __init__(self, **kwargs):

        self.omegab = 0.049
        self.omegac = 0.261
        self.omegam = self.omegab + self.omegac
        self.h      = 0.68
        self.ns     = 0.965
        self.sigma8 = 0.81

        self.c = 3e5

        self.H0 = 100 * self.h

        c  = self.c
        H0 = self.H0
        omegam = self.omegam

        nz = 100000
        z1 = 0.0
        z2 = 6.0

        za = np.linspace(z1,z2,nz)
        dz = za[1]-za[0]

        H      = lambda z: H0*np.sqrt(omegam*(1+z)**3+1-omegam)
        dchidz = lambda z: c/H(z)

        chia = np.cumsum(dchidz(za))*dz

        self.zofchi = interp1d(chia,za)
        self.chiofz = interp1d(za,chia)
