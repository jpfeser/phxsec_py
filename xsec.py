#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 21:39:24 2018


@author: jpfeser
"""
import numpy as np
from numpy import pi,sqrt
from scipy.special import spherical_jn as jsphere, spherical_yn as ysphere

def hsphere(n,z):
    """
    hsphere(n,z) returns the spherical hankel function of the first kind of order n
    defined as jsphere(n,z) + 1j*ysphere(n,z)
    """
    return jsphere(n,z) + 1j*ysphere(n,z)

def GetSigmaSphere(k,p,MatParams):
    """
    Returns:
        sigma: scattering cross section
        scat_eff:  
    Inputs:
        k: a vector of wavenumbers
        p: polarization index (p=1 (longitudinal), p=2,3 (transverse))
        MatParams: Object containing the materials properties
    """
    if p==1:
        sigma = GetSigma_SphereComp(k,MatParams)
    else:
        sigma = GetSigma_SphereTrans(k,MatParams)    
    scat_eff = sigma/(pi*MatParams.a_NP**2)
    return (sigma,scat_eff)

def GetSigma_SphereComp(k,MatParams):
    """
    Calculates scattering cross section due to longitudinal incident waves
    Returns:
        sigma: scattering cross section
    Inputs:
        k: a vector of wavenumbers
        MatParams: Object containing the materials properties
    """
    # First readout some shit from the MatParams input
    a=MatParams.anp
    C11m=MatParams.C11m
    C44m=MatParams.C44m
    rhom=MatParams.rhom
    C11np=MatParams.C11np
    C44np=MatParams.C44np
    rhonp=MatParams.rhonp
    
    vsm = np.array([sqrt(C11m/rhom),sqrt(C44m/rhom),sqrt(C44m/rhom)])
    vsnp = np.array([sqrt(C11np/rhonp),sqrt(C44np/rhonp),sqrt(C44np/rhonp)])
    
    omega = vsm(1)*k
    #k could be a higher dimensional array in principle  
    # So convert it into vector format and we'll reshape the result later.
    omega_vect = omega.flatten() 
    
    omega_zero_logic = (omega_vect==0)
    
    GammaN = TruellXSection_JPFCorrections(omega_vect,a,C11m,C44m,rhom,C11m,C44m,rhom)
    
    sigma = GammaN*pi*a**2 # convert to x-section (m2)
    (n,m)=k.shape # was k a matrix?
    sigma = sigma.reshape(n,m) # make sigma the same shape as k
    return sigma 

def TruellXSection_JPFCorrections(omega_vect,a,C11m,C44m,rhom,C11m,C44m,rhom):
    # still need to be done
    return GammaN

def GetSigma_SphereTrans(k,MatParams):
        """
    Calculates scattering cross section due to transverse incident waves
    Returns:
        sigma: scattering cross section
    Inputs:
        k: a vector of wavenumbers
        MatParams: Object containing the materials properties
    """
    return sigma