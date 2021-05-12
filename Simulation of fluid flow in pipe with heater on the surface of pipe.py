# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 11:54:13 2020

@author: IDEAPAD 310
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


R = 5.1 #cm
Vmax = 0.03 #cm/s
T0 = 30 #C
Nr = 20
dr = R/(Nr-1)
r = np.linspace(0,R,Nr)
L = 60 #cm
a = 1e-3 #cm2/s

T_init = np.ones(Nr)*T0
z = np.linspace(0,L,Nr)


def vz (r):
    return Vmax*(1-(r/R)**2)

def func(T,z):
    dtdz = np.zeros(len(T_init))
    T[-1]=T0+1/20*(60*z-z**2)
    T[0]=(-T[2]+4*T[1])/3
    for i in range (1,Nr-1):
        dtdz [i]=a/vz(r[i])*((T[i+1]-2*T[i]+T[i-1])/dr**2+1/r[i]*(T[i+1]-T[i-1])/(2*dr))
    return dtdz

Tsolv = odeint(func,T_init,z)

Tsolv[:,0]= (-Tsolv[:,2]+4*Tsolv[:,1])/3
for ii in range (len(z)):
    Tsolv[ii][-1]=T0+1/20*(60*z[ii]-z[ii]**2)



plt.figure(0,figsize=(10,10))
plt.imshow(np.transpose(Tsolv),cmap="jet",extent=[0,z[-1],R,0],aspect=z[-1]/R,interpolation='bicubic')
plt.ylabel("radial position")
plt.xlabel("axial position")
plt.colorbar()
