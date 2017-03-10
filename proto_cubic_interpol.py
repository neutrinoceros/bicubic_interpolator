#!/usr/bin/python

#import
import numpy as np
from numpy import random as rd
import matplotlib.pyplot as plt

#defintions
def interpolation_Cubique(p0,p1,p2,p3,x) :
    # dev note ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # I never took into account the x-coordinates of p0,p1,p2,p3
    # this feels wrong on more than one level ...
    # dev note ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    a = -0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3
    b =      p0 - 2.5*p1 + 2.0*p2 - 0.5*p3
    c = -0.5*p0          + 0.5*p2
    d = p1
    return d + c*x + b*x**2 + a*x**3

# def interpolation_BiCubique(p00,p01,p02, ... ,x,y) :
#     f = interpolation_Cubique
#     return f(f(p00,p01,p02,p03,y),
#              f(p10,p11,p12,p13,y),
#              f(p20,p21,p22,p23,y),
#              f(p30,p31,p32,p33,y),x)

#script
N = 4
x = np.linspace(-1,3,N+1)
data = rd.normal(1,1,N+1)

nrad = 1000
refined_x = np.linspace(-1,3,nrad+1)
interpol_data = np.zeros(nrad+1)
i_old = 0
for i in range(nrad) :
    while (x[i_old] < refined_x[i]) :
        i_old += 1
    if i_old > 1 and i_old < N : #forget about the borders for now
        #print i, i_old
        interpol_data[i] = interpolation_Cubique(data[i_old-2],
                                                 data[i_old-1],
                                                 data[i_old],
                                                 data[i_old+1],
                                                 refined_x[i])

plt.scatter(x,data)
#plt.scatter(refined_x,interpol_data, marker='+', color='r')
plt.plot(refined_x,interpol_data, color='r')
plt.savefig("coucou.png")

    
