import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.optimize import curve_fit
import scipy.integrate

from astropy import units as u

import sys
import os

h      = 0.68
omegab = 0.049
omegam = 0.31
omegal = 1-omegam
sigma8 = 0.81
ns     = 0.965

rho_mean    = 2.775e11 * omegam * h**2

# load in power spectrum used for websky
pkfile   = ('./planck2018_powerspectrum.dat')
pk_data = np.loadtxt(pkfile)
k       = pk_data[:,0]
pk      = pk_data[:,1]*(2*np.pi)**3

# routines to get dndm using the Tinker (2008) fitting function
def mass_to_radius(m, rho_mean):
    """
    Returns the radius of a sphere of uniform denstiy rho and mass m

    """
    r = (3 * m / 4. / np.pi / rho_mean)**(1./3.)

    return r

def windowfunction(x):
    """
    Computes the window function in Fourier space (top hat in real space).

    """
    W = (3. / x**3) * (np.sin(x) - x * np.cos(x))

    return W

def M_to_sigma(k, pk, M, omegaM, h):
    """

    Returns mass and sigma(mass)
    """

    rho_mean = 2.775e11 * omegaM * h**2

    # Normalization

    Anorm    = 1./(2.*np.pi**2)

    # Now compute sigma(M)

    sigma = np.zeros(M.shape[0])

    for i in range(sigma.shape[0]):

        radius = mass_to_radius(M[i], rho_mean)

        x = k * radius
        y = pk * (k * windowfunction(x))**2

        sigma[i] = np.sqrt(Anorm * sp.integrate.simps(y, k, even="avg"))

    return sigma

def dlnsigmainv_dM(M, sigma):
    lnsigmainv = np.log(1/sigma)

    diff_sig = np.diff(lnsigmainv)
    diff_M   = np.diff(M)

    dlnsigmainvdM_int = diff_sig / diff_M
    Mint              = M[:-1] + diff_M / 2

    f = interp1d(Mint,dlnsigmainvdM_int,fill_value="extrapolate")

    return f(M)

def growth_factor(omegaM, omegaL, z):
    """
    Returns growth factor using fitting formulae from Carrol, Press & Turner (1992)
    """
    # returns growth factor using fitting formulae from Carrol, Press & Turner (1992)

    w = -1.
    x = 1.+z
    x2 = x**2
    x3 = x**3
    x3w = x**(3*w)

    #calculate omegaM(z) and omegaL(z)
    denom = omegaM*x3 + (1-omegaM-omegaL)*x2 + omegaL*x3*x3w
    omega = omegaM*x3/denom
    lamb = omegaL*x3*x3w/denom

    #fitting formulae for w=-1 universe
    g = 2.5*omega/(omega**(4.0/7.0) - lamb + (1+(omega/2))*(1+(lamb/70)))
    g0 = 2.5*omegaM/(omegaM**(4.0/7.0) - omegaL + (1+(omegaM/2))*(1+(omegaL/70)))

    D = (g/x)/g0

    return D

def tinker_func(x, z):
    """
    Uses fitting coefficients from Table 2 of Tinker et al. (2008) for Delta = 200 and
    redshift evolution equations 5 through 8.
    """

    A_hmfcalc = 1.858659e-01
    a_hmfcalc = 1.466904
    b_hmfcalc = 2.571104
    c_hmfcalc = 1.193958

    A =  0.186
    a = 1.47
    b = 2.57
    c = 1.19

    A =  A_hmfcalc
    a =  a_hmfcalc
    b =  b_hmfcalc
    c =  c_hmfcalc

    amp = A * (1. + z)**-0.14
    a = a * (1. + z)**-0.06
    alpha = 10**(-1. * (0.75 / np.log10(200. / 75.))**1.2)

    b = b * (1. + z)**-alpha
    c = c

    f = amp * ((x / b)**(-a) + 1.) * np.exp(-c / x**2)

    return f

def dndmofm_tinker(Mmin, Mmax, redshift):
    """
    Returns dn/dm for Tinker et al. (2008) mass function using Eqs. (3, 5-8) of their paper.
    Assumes k, pk, omegam, omegal, h, and rho_mean have been defined already
    units: m [Msun], dndm [1/Msun/Mpc^3]
    """

    dlog = 0.001
    n    = int((np.log10(Mmax)-np.log10(Mmin))/dlog)
    M    = np.logspace(np.log10(Mmin),np.log10(Mmax),n)

    sigma  = M_to_sigma(k, pk, M, omegam, h)
    D      = growth_factor(omegam, omegal, redshift)
    sigma *= D

    fsigma = tinker_func(sigma, redshift)

    dlnsigmainv = dlnsigmainv_dM(M, sigma)

    dndm = tinker_func(sigma, redshift) * rho_mean / M * dlnsigmainv

    dndmofm = interp1d(M,dndm,fill_value='extrapolate')

    return M,dndm,dndmofm

def dndmofmz_tinker(mmin,mmax,zmin,zmax):

    from os import path
    if not path.exists('dndmtab.npz'):

        print("\n creating table")
        dlogm = 0.05
        dz    = 0.1

        nm = int((np.log10(mmax)-np.log10(mmin))/dlogm)+1
        nz = int((zmax-zmin)/dz)+1

        m = np.logspace(np.log10(mmin),np.log10(mmax),nm)
        z = np.linspace(zmin,zmax,nz)

        dndmofmz = np.zeros((nz,nm))

        iz=0
        for zv in z:
            print(" z: ","{:4.2f}".format(zv),
                  end="\r", flush=True)
            m_t,dndm_t, f = dndmofm_tinker(mmin,mmax,zv) # M, dndM in Msun, 1/Mpc^3/Msun
            dndmofmz[iz,:] = np.log10(f(m))
            iz += 1
        np.savez('dndmtab.npz',m=m,z=z,dndmofmz=dndmofmz)
    else:
        data = np.load('dndmtab.npz')
        m = data['m']
        z = data['z']
        dndmofmz = data['dndmofmz']

    dndmofmzfunc_log = interp2d(np.log10(m),z,dndmofmz)
    dndmofmzfunc     = lambda m,z: 10**dndmofmzfunc_log(np.log10(m),z)

    return dndmofmzfunc
