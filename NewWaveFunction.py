#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 14:47:22 2019

@author: physics_001
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, sqrt, atan
pi = np.pi
import statistics as stats
from scipy.signal import find_peaks_cwt


iter = 140
x = np.arange(-10,11,1) #x-array
y = np.arange(-10,11,1) #y-array
notlist = np.arange(1,iter + 1,1)
thetaA = (2*pi/iter)*(notlist)# list of angles that is the same size as x and y array
results = []

#WaveFunction
def WF(x,y):
    return cos(pi/4)*x + cos(pi/4)*y
#Gaussian
def Gaussian(x, y, theta, x0 = 0, y0 = 0, sigma =.1):
    return np.exp(-(((x - x0)**2 + (y - y0)**2))/4*sigma**2 + (1j)*(cos(theta)*x + sin(theta)*y))

# make an empty matrix of the size of the x and y vectors for the WaveFunction and the Gaussian
F = np.zeros([len(y),len(x)], dtype=complex)
G = np.zeros([len(y),len(x)], dtype=complex)

i = 0
j = 0
count = 0
for k in thetaA:
    for j in range(len(y)):
        for i in range(len(x)):
            theta = thetaA[count]
            F[j,i] = WF(x[i], y[j])
            G[j,i] = Gaussian(x[i], y[j], theta)
    Z = np.sum(F*G, axis = 1)
    D = np.sum(Z)
    mag = sqrt((D.real)**2 + (D.imag)**2)
    results.append(float(mag))
    count = count + 1

#print(results)
maxr = max(results)
r = []
for value in range(len(results)):
    r.append(results[value]/maxr)
    
    
def compass(x, y, arrowprops=None):
    angles = thetaA
    radii = r
    fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(15,15),)
    
    ax.grid(False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.spines['polar'].set_visible(False)
    kw = dict(arrowstyle="->", color='black')
    if arrowprops:
        kw.update(arrowprops)
    [ax.annotate("", xy=(angle, radius), xytext=(0, 0),
                 arrowprops=kw) for
     angle, radius in zip(angles, radii)]
   

    ax.set_ylim(0, np.max(radii))

    return fig, ax

fig, ax = compass(theta, r)

plt.plot(thetaA,results)
plt.show()
