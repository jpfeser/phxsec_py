#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 17:29:05 2018

@author: jpfeser
"""
import numpy as np
from matplotlib import pyplot as plt
import xsec

anp = 3e-9
C11m = 105.7e9
C44m = 44e9
rhom = 5500
C11np = 158.5e9
C44np = 40e9
rhonp = 8250

#  Ge in Si
anp = 3e-9

C11m = 166e9
C44m = 80e9
rhom = 2300

C11np = 126e9
C44np = 68e9
rhonp = 5323

MP = xsec.create_matparams(anp,C11m,C44m,rhom,C11np,C44np,rhonp)

ka_vect = np.logspace(-1.0,2.0,num=200)
#(a,C11m,C44m,rhom,C11np,C44np,rhonp)=xsec.readout_matparams(MatParams)
(sigma, scat_eff) = xsec.sigma_sphere(ka_vect/anp,2,MP)
(sigma_long, scat_eff_long) = xsec.sigma_sphere(ka_vect/anp,1,MP)

plt.figure(1)
plt.subplot(211)
plt.loglog(ka_vect,scat_eff, label = 'transverse')
plt.loglog(ka_vect,scat_eff_long, label = 'longitudinal')
#plt.savefig('output.eps',format='eps')
plt.xlabel(r'$ka$')
plt.ylabel(r'$\chi$')
plt.legend()
plt.show()

plt.subplot(212)
plt.loglog(ka_vect,4/3*ka_vect/0.5*1/scat_eff, label = 'transverse')
plt.loglog(ka_vect,4/3*ka_vect/0.5*1/scat_eff_long, label = 'longitudinal')
#plt.savefig('output.eps',format='eps')
plt.xlabel(r'$ka$')
plt.ylabel(r'$k\ell$')
plt.legend()
plt.show()