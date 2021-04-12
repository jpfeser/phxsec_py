#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 21:39:24 2018


@author: jpfeser
"""
import numpy as np
from numpy import pi,sqrt
from scipy.special import spherical_jn as jsphere, spherical_yn as ysphere

# def __init__:
class MatParams:
    anp = 3e-9
    C11m = 80e9
    C44m = 40e9
    rhom = 5500
    C11np = 300e9
    C44np = 40e9
    rhonp = 5600
    
def hsphere(n,z):
    """
    hsphere(n,z) returns the spherical hankel function of the first kind of order n
    defined as jsphere(n,z) + 1j*ysphere(n,z)
    """
    return jsphere(n,z) - 1j*ysphere(n,z)

def readout_matparams(MatParams):
    """
    Returns:
        the tuple: (anp,C11m,C44m,rhom,C11np,C44np,rhonp)  
    Inputs:
        MatParams: Object containing the materials properties MatParms.anp etc
    Example:  (anp,C11m,C44m,rhom,C11np,C44np,rhonp)=readout_matparams(MatParams)
    """
    anp=MatParams.anp
    C11m=MatParams.C11m
    C44m=MatParams.C44m
    rhom=MatParams.rhom
    C11np=MatParams.C11np
    C44np=MatParams.C44np
    rhonp=MatParams.rhonp
    
    return (anp,C11m,C44m,rhom,C11np,C44np,rhonp)

    """
    Returns:
        the tuple: (anp,C11m,C44m,rhom,C11np,C44np,rhonp)  
    Inputs:
        MatParams: Object containing the materials properties MatParms.anp etc
    Example:  (anp,C11m,C44m,rhom,C11np,C44np,rhonp)=readout_matparams(MatParams)
    """
class create_matparams:
    def __init__(self, anp,C11m,C44m,rhom,C11np,C44np,rhonp):
        self.anp = anp
        self.C11m = C11m
        self.C44m = C44m
        self.rhom = rhom
        self.C11np = C11np
        self.C44np = C44np
        self.rhonp = rhonp

def sigma_sphere(k,p,MatParams):
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
        sigma = sigma_sphere_comp(k,MatParams)
    else:
        sigma = sigma_sphere_trans(k,MatParams)
    (a,C11m,C44m,rhom,C11np,C44np,rhonp)=readout_matparams(MatParams)
    scat_eff = sigma/(pi*a**2)
    return (sigma,scat_eff)

def sigma_sphere_comp(k,MatParams):
    """
    Calculates scattering cross section due to longitudinal incident waves
    Returns:
        sigma: scattering cross section
    Inputs:
        k: a vector of wavenumbers
        MatParams: Object containing the materials properties
    """
    # readout Materials Properties
    (a,C11m,C44m,rhom,C11np,C44np,rhonp)=readout_matparams(MatParams)
    
    # construct sound speed vector: [v_long,v_trans,v_trans] in matrix
    vsm = np.array([sqrt(C11m/rhom),sqrt(C44m/rhom),sqrt(C44m/rhom)])
    
    # calculate omega array (rad/s) from k array
    omega_in = vsm[0]*k
    #k could be a higher dimensional array in principle  
    # So convert it into vector format and we'll reshape the result later.
    omega = omega_in.flatten() 
    
    #
    mmax = 200
    nomega = omega.size
    scat_eff= np.zeros_like(omega)
    for n in np.arange(nomega):
        k1 = omega[n]*(C11m/rhom)**(-0.5)
        kap1 = omega[n]*(C44m/rhom)**(-0.5)
        k2 = omega[n]*(C11np/rhonp)**(-0.5)
        kap2 = omega[n]*(C44np/rhonp)**(-0.5)
        
        k1a = k1*a
        kap1a = kap1*a
        k2a = k2*a
        kap2a = kap2*a
        
        kapokap = kap1/kap2
        ror = rhom/rhonp
        mom = C44m/C44np
        
        #First calculate the m=0 coeficients
        Eterm = 1/ror*((kap1a)**2/(1-k2a*np.tan(k2a)**(-1))-4*(kapokap)**2)
        FE= 4+Eterm
        num = ((kap1a)**2-FE)*np.sin(k1a)+FE*k1a*np.cos(k1a)
        denom = ((kap1a**2-FE)**2+FE**2*k1a**2)**0.5
        multiplier = np.exp(1j*(k1a-np.arctan((kap1a**2-FE)/(FE*k1a))))
        
        A=np.zeros((mmax,1),dtype='complex64')
        B=np.zeros((mmax,1),dtype='complex64')
        C=np.zeros((mmax,1),dtype='complex64')
        D=np.zeros((mmax,1),dtype='complex64')
        A0_Truell = 1/(k1a)*num/denom*multiplier
    
    #First calculate the m=0 coeficients
        m=0
        AA = np.zeros((8,8))
        delta = np.zeros((8,1))
   
        AA[0,0] = k1a*jsphere(m+1,k1a)
        AA[0,1] = k1a*ysphere(m+1,k1a)
        AA[0,2] = 0
        AA[0,3] = 0
        AA[0,4] = -k2a*jsphere(m+1,k2a)
        AA[0,5] = 0
        AA[0,6] = 0
        AA[0,7] = 0
        
        AA[2,0] = jsphere(m,k1a)
        AA[2,1] = ysphere(m,k1a)
        AA[2,2] = -(jsphere(m,kap1a)-kap1a*jsphere(m+1,kap1a))
        AA[2,3] = -(ysphere(m,kap1a)-kap1a*ysphere(m+1,kap1a))
        AA[2,4] = -jsphere(m,k2a)
        AA[2,5] = 0
        AA[2,6] = jsphere(m,kap2a)-kap2a*jsphere(m+1,kap2a)
        AA[2,7] = 0
        
        AA[4,0] = kap1a**2*jsphere(m,k1a)-2*(m+2)*k1a*jsphere(m+1,k1a)
        AA[4,1] = kap1a**2*ysphere(m,k1a)-2*(m+2)*k1a*ysphere(m+1,k1a)
        AA[4,2] = 0
        AA[4,3] = 0
        AA[4,4] = -1/mom*(kap2a**2*jsphere(m,k2a)-2*(m+2)*k2a*jsphere(m+1,k2a))
        AA[4,5] = 0
        AA[4,6] = 0
        AA[4,7] = 0
        
        AA[6,0] = (m-1)*jsphere(m,k1a)-k1a*jsphere(m+1,k1a)
        AA[6,1] = (m-1)*ysphere(m,k1a)-k1a*ysphere(m+1,k1a)
        AA[6,2] = -((m**2-1-kap1a**2/2)*jsphere(m,kap1a)+kap1a*jsphere(m+1,kap1a))
        AA[6,3] = -((m**2-1-kap1a**2/2)*ysphere(m,kap1a)+kap1a*ysphere(m+1,kap1a))
        AA[6,4] = -1/mom*((m-1)*jsphere(m,k2a)-k2a*jsphere(m+1,k2a))
        AA[6,5] = 0
        AA[6,6] = 1/mom*((m**2-1-kap2a**2/2)*jsphere(m,kap2a)+kap2a*jsphere(m+1,kap2a))
        AA[6,7] = 0
        
        for kk in np.arange(1,5):  #1 through 4
            for LL in np.arange(1,9): #1 through 8
                #odd number
                if np.mod(LL,2)==1:
                    AA[2*kk-1,LL-1] = -AA[2*kk-2,LL]
                #even numbers
                else:  
                    AA[2*kk-1,LL-1] = AA[2*kk-2,LL-2]
                    
        for LL in np.arange(1,9): #1 through 8
            #odd number
            if np.mod(LL,2)==1:
                delta[LL-1] = (-1)**m*(1/k1a)*AA[LL-1,0]
            #even numbers
        else:
            delta[LL-1] = 0
            
        #Now solve
        result = np.linalg.lstsq(AA,delta,rcond=None)
        X = result[0]
        residuals = result[1]
        rank = result[2]
        s = result[3]
        
        A[m] = X[0] + 1j*X[1]
        B[m] = X[2] + 1j*X[3]
        C[m] = X[4] + 1j*X[5]
        D[m] = X[6] + 1j*X[7]
        
        summand = np.zeros((mmax,1))
        summand[0] = 4*(np.abs(A[0]))**2
        
        #awsome, m=0 term complete.
        reltol=1e-10
        relerr=1
        while relerr>reltol:
            m=m+1
            AA = np.zeros((8,8))
            delta = np.zeros((8,1))
       
            AA[0,0] = k1a*jsphere(m+1,k1a)-m*jsphere(m,k1a)
            AA[0,1] = k1a*ysphere(m+1,k1a)-m*ysphere(m,k1a)
            AA[0,2] = m*(m+1)*jsphere(m,kap1a)
            AA[0,3] = m*(m+1)*ysphere(m,kap1a)
            AA[0,4] = -k2a*jsphere(m+1,k2a)+m*jsphere(m,k2a)
            AA[0,5] = 0
            AA[0,6] = -m*(m+1)*jsphere(m,kap2a)
            AA[0,7] = 0
            
            AA[2,0] = jsphere(m,k1a)
            AA[2,1] = ysphere(m,k1a)
            AA[2,2] = -((m+1)*jsphere(m,kap1a)-kap1a*jsphere(m+1,kap1a))
            AA[2,3] = -((m+1)*ysphere(m,kap1a)-kap1a*ysphere(m+1,kap1a))
            AA[2,4] = -jsphere(m,k2a)
            AA[2,5] = 0
            AA[2,6] = (m+1)*jsphere(m,kap2a)-kap2a*jsphere(m+1,kap2a)
            AA[2,7] = 0
            
            AA[4,0] = kap1a**2*jsphere(m,k1a)-2*(m+2)*k1a*jsphere(m+1,k1a)
            AA[4,1] = kap1a**2*ysphere(m,k1a)-2*(m+2)*k1a*ysphere(m+1,k1a)
            AA[4,2] = m*(kap1a**2*jsphere(m,kap1a)-2*(m+2)*kap1a*jsphere(m+1,kap1a))
            AA[4,3] = m*(kap1a**2*ysphere(m,kap1a)-2*(m+2)*kap1a*ysphere(m+1,kap1a))
            AA[4,4] = -1/mom*(kap2a**2*jsphere(m,k2a)-2*(m+2)*k2a*jsphere(m+1,k2a))
            AA[4,5] = 0
            AA[4,6] = -1/mom*m*(kap2a**2*jsphere(m,kap2a)-2*(m+2)*kap2a*jsphere(m+1,kap2a))
            AA[4,7] = 0
            
            AA[6,0] = (m-1)*jsphere(m,k1a)-k1a*jsphere(m+1,k1a)
            AA[6,1] = (m-1)*ysphere(m,k1a)-k1a*ysphere(m+1,k1a)
            AA[6,2] = -((m**2-1-kap1a**2/2)*jsphere(m,kap1a)+kap1a*jsphere(m+1,kap1a))
            AA[6,3] = -((m**2-1-kap1a**2/2)*ysphere(m,kap1a)+kap1a*ysphere(m+1,kap1a))
            AA[6,4] = -1/mom*((m-1)*jsphere(m,k2a)-k2a*jsphere(m+1,k2a))
            AA[6,5] = 0
            AA[6,6] = 1/mom*((m**2-1-kap2a**2/2)*jsphere(m,kap2a)+kap2a*jsphere(m+1,kap2a))
            AA[6,7] = 0
            
            for kk in np.arange(1,5):  #1 through 4
                for LL in np.arange(1,9): #1 through 8
                    #odd number
                    if np.mod(LL,2)==1:
                        AA[2*kk-1,LL-1] = -AA[2*kk-2,LL]
                    #even numbers
                    else:  
                        AA[2*kk-1,LL-1] = AA[2*kk-2,LL-2]
                        
            for LL in np.arange(1,9): #1 through 8
                #odd number
                if np.mod(LL,2)==1:
                    delta[LL-1] = (-1)**m*(1/k1a)*AA[LL-1,0]
                #even numbers
            else:
                delta[LL-1] = 0
                
            #Now solve
            result = np.linalg.lstsq(AA,delta,rcond=None)
            X = result[0]
            residuals = result[1]
            rank = result[2]
            s = result[3]
            
            A[m] = X[0] + 1j*X[1]
            B[m] = X[2] + 1j*X[3]
            C[m] = X[4] + 1j*X[5]
            D[m] = X[6] + 1j*X[7]
            
            summand[m] = 4*(2*m+1)*((np.abs(A[m]))**2+m*(m+1)*(k1/kap1)*(np.abs(B[m]))**2)
            relerr = np.abs(summand[m])/np.sum(summand)
            #print(relerr)
            
        scat_eff[n]=np.sum(summand)
            
    sigma = scat_eff*pi*a**2 # convert to x-section (m2)
    if np.isscalar(k):
        return sigma 
    else:
        # k was an array.  
        #Make sure things go out the same form they came in. 
        sigma = sigma.reshape(k.shape)
        return sigma 

def sigma_sphere_trans(k,MatParams):
    """
    Calculates scattering cross section due to transverse incident waves
    Returns:
        sigma: scattering cross section
    Inputs:
        k: a vector of wavenumbers
        MatParams: Object containing the materials properties
    """
    (a,C11_1,C44_1,rho_1,C11_2,C44_2,rho_2)=readout_matparams(MatParams)
    
#nomenclature/calculation is based on:
#Iwashimazu, Journal of Sound and Vibration, 40(2), 267-271, 1975
#
# NOTE: subscript 2 is the Nanoparticle in this code!!!

    # construct sound speed vector: [v_long,v_trans,v_trans] in matrix
    cL_1 = sqrt(C11_1/rho_1)
    cT_1 = sqrt(C44_1/rho_1)
    cL_2 = sqrt(C11_2/rho_2)
    cT_2 = sqrt(C44_2/rho_2)
    
    # Lame constants
    mu_1 = C44_1
    lambda_1 = C11_1 - 2*mu_1
    mu_2 = C44_2
    lambda_2 = C11_2 - 2*mu_2
    
    # calculate omega array (rad/s) from k array
    omega_in = cT_1*k
    #k could be a higher dimensional array in principle  
    # So convert it into vector format and we'll reshape the result later.
    omega = omega_in.flatten() 
    
    #
    mmax = 200
    QL = np.zeros_like(omega)
    QT1 =np.zeros_like(omega)
    QT2 =np.zeros_like(omega)
    nomega = omega.size
    scat_eff= np.zeros_like(omega)
    for no in np.arange(nomega):
        # print(no+1)
        Ka_L_1 = omega[no]/cL_1*a
        Ka_T_1 = omega[no]/cT_1*a
    
        Ka_L_2 = omega[no]/cL_2*a
        Ka_T_2 = omega[no]/cT_2*a
        
        p = mu_2/mu_1
        kappa = Ka_L_1/Ka_T_1
        
        QTot_loop = 0
        QLn = np.zeros((mmax,1))
        QT1n = np.zeros((mmax,1))
        QT2n = np.zeros((mmax,1))
        n=0
        reltol=1e-10
        relerr=1
        while relerr>reltol:
            n = n + 1
            J_ka_1 = Ka_T_1*jsphere(n+1,Ka_T_1)/jsphere(n,Ka_T_1)
            H_ka_1 = Ka_T_1*hsphere(n-1,Ka_T_1)/hsphere(n,Ka_T_1)
            J_Ka_1 = Ka_L_1*jsphere(n+1,Ka_L_1)/jsphere(n,Ka_L_1)
            H_Ka_1 = Ka_L_1*hsphere(n-1,Ka_L_1)/hsphere(n,Ka_L_1)
            J_ka_2 = Ka_T_2*jsphere(n+1,Ka_T_2)/jsphere(n,Ka_T_2)
            H_ka_2 = Ka_T_2*hsphere(n-1,Ka_T_2)/hsphere(n,Ka_T_2)
            J_Ka_2 = Ka_L_2*jsphere(n+1,Ka_L_2)/jsphere(n,Ka_L_2)
            H_Ka_2 = Ka_L_2*hsphere(n-1,Ka_L_2)/hsphere(n,Ka_L_2)
            
            #column vector
            A1 = np.array([[-(n+1) + H_Ka_1],
                            [1],
                            [2*(n+1)*(n+2)-((Ka_T_1)**2+4*H_Ka_1)],
                            [-2*(n+2)+2*H_Ka_1]])
            #column vector
            A2 = np.array([[-(n+1)],
                    [ 1-H_ka_1/n ],
                    [ 2*(n+1)*(n+2)-2*(n+1)*H_ka_1 ],
                    [ -2*(n+2)+((Ka_T_1)**2+2*H_ka_1)/n]])
    
            #column vector
            A3 = np.array([[n-J_Ka_1],
                           [1],
                           [2*n*(n-1)*p + (-Ka_T_2**2 + 4*J_Ka_2)*p],
                           [2*(n-1)*p - 2*J_Ka_2*p]])
            #column vector
            A4 = np.array([[n],
                [1-J_ka_2/(n+1)],
                [2*n*(n-1)*p-2*n*J_ka_2*p],
                [2*(n-1)*p+(-Ka_T_2**2 + 2*J_ka_2)*p/(n+1)]])
    
            B = np.array([[1],
                          [1/n - J_ka_1/(n*(n+1))],
                          [2*(n-1)-2*J_ka_1],
                          [2*(n-1)/n-((Ka_T_1)**2-2*J_ka_1)/(n*(n+1))]])
    
            C1 = np.array([[1],
                  [-(n+2) + H_ka_1]])
    
            C2 = np.array([[1],
                           [(n-1)*p - p*J_ka_2]])
    
            D = np.array([[1],
                          [(n-1)-J_ka_1]])
    
            detA = np.linalg.det(np.block([A1,A2,A3,A4]))
            detC = np.linalg.det(np.block([C1,C2]))
            det1 = np.linalg.det(np.block([B,A2,A3,A4]))
            det2 = np.linalg.det(np.block([A1,B,A3,A4]))
            det3 = np.linalg.det(np.block([D,C2]))
            det4 = np.linalg.det(np.block([A1,A2,B,A4]))
            det5 = np.linalg.det(np.block([A1,A2,A3,B]))
            det6 = np.linalg.det(np.block([C1,D]))
            
            phi_1 = 1j**(n+1)*(2*n+1)*(jsphere(n,Ka_T_1)/(Ka_T_1*hsphere(n,Ka_L_1)))*det1/detA
            psi_1 = 1j**(n+1)*(2*n+1)/n*(jsphere(n,Ka_T_1)/(Ka_T_1*hsphere(n,Ka_T_1)))*det2/detA
            omega_1 = 1j**(n)*(2*n+1)/(n*(n+1))*(jsphere(n,Ka_T_1)/(hsphere(n,Ka_T_1)))*det3/detC
            
            phi_2 = -1j**(n+1)*(2*n+1)*(jsphere(n,Ka_T_1)/(Ka_T_1*jsphere(n,Ka_L_2)))*det4/detA
            psi_2 = 1j**(n+1)*((2*n+1)/n+1)*(jsphere(n,Ka_T_1)/(Ka_T_1*jsphere(n,Ka_T_2)))*det5/detA
            omega_2 = -1j**(n)*(2*n+1)/(n*(n+1))*(jsphere(n,Ka_T_1)/(jsphere(n,Ka_T_2)))*det6/detC;
            
            pref = 2*n*(n+1)/(2*n+1)
            QLn[n-1] = pref*(Ka_T_1/Ka_L_1)*np.abs(phi_1)**2 
            QT1n[n-1] = pref*n*(n+1)*np.abs(psi_1)**2
            QT2n[n-1] = pref*n*(n+1)*np.abs(omega_1)**2/(Ka_T_1)**2
        
            QTotprev = QTot_loop
            QTot_loop = QTot_loop + (QT1n[n-1] + QT2n[n-1] + QLn[n-1])
            if n>1:
                relerr = (QTot_loop-QTotprev)/QTotprev
            # end while loop
            
        QT1[no]=np.sum(np.real(QT1n))
        QT2[no]=np.sum(np.real(QT2n))
        QL[no]=np.sum(np.real(QLn))
        #end for loop over different omega values
        
    scat_eff = QT1 + QT2 + QL
    sigma = scat_eff*pi*a**2
    if np.isscalar(k):
        return sigma 
    else:
        # k was an array.  
        #Make sure things go out the same form they came in. 
        sigma = sigma.reshape(k.shape)
        return sigma 