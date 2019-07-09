#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 08:43:17 2019

@author: physics_001
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, atan, sqrt
pi = np.pi
import statistics as stat
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy import integrate

###############################################################################
## Generate Data for Wavefunction and Gaussian.

iter = 240
x = np.arange(-10,11,1) #x-array
y = np.arange(-10,11,1) #y-array
increment = np.arange(1,iter+1,1)

thetaA = (2*pi/iter)*(increment) # Doesn't include 0 rad?
results = []

def WF(x,y):
    return (1*cos(pi/4)*x + 1*sin(pi/4)*y)

#def WF(x,y):
#    return np.exp(1j*(cos(pi/4)*x + sin(pi/4)*y))

def Gaussian(x, y, theta, x0 = 0, y0 = 0, sigma = 0.1):
    return np.exp(-(((x - x0)**2 + (y - y0)**2))/4*sigma**2 + (1j)*(cos(theta)*x + sin(theta)*y))

# make an empty matrix of the size of the x and y vectors for the WaveFunction and the Gaussian
F = np.zeros([len(y),len(x)], dtype=complex)
G = np.zeros([len(y),len(x)], dtype=complex)

i = 0
j = 0

count = 0


for k in thetaA:
    #print(k)
    #print(count)
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
#print(stat.mode(results))
#plt.plot(results)
#print(results)
maxr = max(results)
r = []
for value in range(len(results)):
    r.append(results[value]/maxr)

###############################################################################
## Find local Max and Min for r.
  
a = np.diff(np.sign(np.diff(r))).nonzero()[0] + 1            # local min & max
mins = (np.diff(np.sign(np.diff(r))) > 0).nonzero()[0] + 1         # local min
maxs = (np.diff(np.sign(np.diff(r))) < 0).nonzero()[0] + 1         # local max
# +1 due to the fact that diff reduces the original index number

###############################################################################
## "Shift" data to start at a minimum.

thetaAadd = (2*pi/iter)*(increment) # Doesn't include 0 rad?
for angle in range(mins[0]):
    thetaAadd[angle] = thetaAadd[angle] + (2*pi)
thetaAnew = np.roll(thetaAadd, (thetaAadd.size - mins[0]))
rnew = r + r[0:mins[0]] # doesn't include value associated with mins[0] i.e. first value is now first minimum.
del rnew[0:mins[0]] # deletes the data from the beginning that was added to the end.
#print(thetaA[mins[0]], r[mins[0]])
###############################################################################
## Find local Max and min for rnew

#a = np.diff(np.sign(np.diff(r))).nonzero()[0] + 1            # local min & max
minsn = []
for newmin in range(len(mins)):
    minsn.append(mins[newmin] - mins[0])        # local min
maxsn = []
for newmax in range(len(maxs)):
    maxsn.append(maxs[newmax] - maxs[0])         # local max
# +1 due to the fact that diff reduces the original index number


###############################################################################

    
plt.figure(figsize=(10, 5))
plt.plot(thetaAnew, rnew, color='grey')
#plt.hist(rnew, bins=240, range=[thetaAnew[0], thetaAnew[len(thetaAnew)-1]])
#plt.plot(thetaA, r, color='green')
for bat in minsn:
    plt.plot(thetaAnew[bat], rnew[bat], "o", label="minimums", color='r')
for cat in maxs:
    plt.plot(thetaA[cat], r[cat], "o", label="maximums", color='b')
#plt.show() 
  
###############################################################################
## Calculate processed Husimi vectors. (Average magnitude and angle)
r_avg = []
theta_avg = []
left = 0
right = left + 1
for peak in range(len(minsn)): # changed from maxsn
#    print(left, right)
    if left == len(minsn)-1:
        ## Ensure that the last value in each dataset is used for the final range.
        ravg = stat.mean(rnew[minsn[left]:-1])
        r_avg.append(ravg)
        thetaavg = stat.mean(thetaAnew[minsn[left]:-1])
        theta_avg.append(thetaavg)
#        plt.plot(thetaA[mins[left]:-1], rnew[mins[left]:-1] )
    else:
        ravg = stat.mean(rnew[minsn[left]:minsn[right]])
        r_avg.append(ravg)
        thetaavg = stat.mean(thetaAnew[minsn[left]:minsn[right]])
        theta_avg.append(thetaavg)
#        plt.plot(thetaA[mins[left]:mins[right]], rnew[mins[left]:mins[right]])
##     write code here to calc. Husimi vector.
        left = left + 1
        right = right + 1
#print(len(r_avg))
for value in range(len(r_avg)):
    plt.plot(theta_avg[value], r_avg[value], "o", label="local avg.", color='orange')
#plt.legend(loc='upper left', ncol=2, shadow=False)
plt.show()
#print(minsn)

#print(thetaAnew[mins[1]], )
###############################################################################
## Could add additional if statement that skipped small ranges of minimums
##      to include other peaks with the larger more defined peaks.        

## Once the processed Husimi (proc-H) map is obtained above, the following lines will find 
##      the peaks of the proc-H map and produce a second proc-H map. 
  
## Threshold can be set after taking the average to remove any Husimi vectors
##      contributing a given amount.


###############################################################################

## Plot Husimi map and proc-H map.    
"""
def compass(u, v, arrowprops=None):
    angles = thetaA
    radii = r
    fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(8,8))
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.spines['polar'].set_visible(False)
    kw = dict(arrowstyle="->", color='gray')
    if arrowprops:
        kw.update(arrowprops)
    [ax.annotate("", xy=(angle, radius), xytext=(0, 0),
                 arrowprops=kw) for
     angle, radius in zip(angles, radii)]
    kw2 = dict(arrowstyle="->", color='red')
    [ax.annotate("", xy=(angle, radius), xytext=(0, 0),
                 arrowprops=kw2) for
     angle, radius in zip(theta_avg, r_avg)]  
    ax.set_ylim(0, np.max(radii))

    return fig, ax

fig, ax = compass(theta, r)

#for vector in range(len(r_avg)):    
#    plt.arrow(theta_avg[vector], 0, 0, r_avg[vector], width = 0.001,
#              edgecolor = 'red', facecolor = 'red', head_width=0.055)   
    
###############################################################################    
  

















