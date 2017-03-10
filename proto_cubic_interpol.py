#!/usr/bin/python

#import
import numpy as np
from numpy import random as rd
import matplotlib.pyplot as plt

# defintions -------------------------------------------------------------
def update_coefficients(x0,p0,x1,p1,x2,p2,x3,p3) :
    # those formulas where calculated and copied/pasted from Mathematica
    a=p0/((x0-x1)*(x0-x2)*(x0-x3))\
       -p1/((x0-x1)*(x1-x2)*(x1-x3))\
       +p2/((x0-x2)*(x1-x2)*(x2-x3))\
       -p3/((x0-x3)*(-x1+x3)*(-x2+x3))
    b=(p2*(x0+x1+x3))/((x0-x2)*(-x1+x2)*(x2-x3))\
        +(p3*(x0+x1+x2))/((x0-x3)*(-x1+x3)*(-x2+x3))\
        +(p1*(x0+x2+x3))/((x0-x1)*(x1-x2)*(x1-x3))\
        -(p0*(x1+x2+x3))/((x0-x1)*(x0-x2)*(x0-x3))
    c=-((p3*(x1*x2+x0*(x1+x2)))/((x0-x3)*(-x1+x3)*(-x2+x3)))-\
        (p2*(x1*x3+x0*(x1+x3)))/((x0-x2)*(-x1+x2)*(x2-x3))\
        -(p1*(x2*x3+x0*(x2+x3)))/((x0-x1)*(x1-x2)*(x1-x3))\
        +(p0*(x2*x3+x1*(x2+x3)))/((x0-x1)*(x0-x2)*(x0-x3))
    d=-((p0*x1*x2*x3)/((x0-x1)*(x0-x2)*(x0-x3)))\
        +x0*((p1*x2*x3)/((x0-x1)*(x1-x2)*(x1-x3))\
             +(p2*x1*x3)/((x0-x2)*(-x1+x2)*(x2-x3))\
             +(p3*x1*x2)/((x0-x3)*(-x1+x3)*(-x2+x3)))
    return a,b,c,d

def interpolation_Cubique(a,b,c,d,x) :
    return d + c*x + b*x**2 + a*x**3

# def interpolation_BiCubique(p00,p01,p02, ... ,x,y) :
#     f = interpolation_Cubique
#     return f(f(p00,p01,p02,p03,y),
#              f(p10,p11,p12,p13,y),
#              f(p20,p21,p22,p23,y),
#              f(p30,p31,p32,p33,y),x)


# script -----------------------------------------------------------------
N = 10
x = np.linspace(-1,2,N+1)
data = rd.normal(1,0.1,N+1)

nrad = 1000
refined_x = np.linspace(-1,2,nrad+1)
interpol_data = np.zeros(nrad+1)

A,B,C,D = 0.,0.,0.,0.
i_old = 0
# uncomment this when ready to dev border conditions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#X0,X1,X2,X3 = x   [i_old-2:i_old+2]
#P0,P1,P2,P3 = data[i_old-2:i_old+2]
#A,B,C,D     = update_coefficients(X0,P0,X1,P1,X2,P2,X3,P3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for i in range(nrad) :
    update_required = False
    while (x[i_old] < refined_x[i]) :
        i_old += 1
        update_required = True
    if update_required and i_old > 1 and i_old < N:
        X0,X1,X2,X3 = x   [i_old-2:i_old+2]
        P0,P1,P2,P3 = data[i_old-2:i_old+2]
        A,B,C,D = update_coefficients(X0,P0,X1,P1,X2,P2,X3,P3)
    if i_old > 1 and i_old < N : #forget about the borders for now
        interpol_data[i] = interpolation_Cubique(A,B,C,D,refined_x[i])

plt.scatter(x,data)
plt.plot(refined_x,interpol_data, color='r')
plt.savefig("coucou.png")

    
