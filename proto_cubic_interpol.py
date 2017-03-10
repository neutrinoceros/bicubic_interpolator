#!/usr/bin/python

# imports ----------------------------------------------------------------

import numpy as np
from numpy import random as rd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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
    return a*x**3 + b*x**2 + c*x + d

# def interpolation_BiCubique(p00,p01,p02, ... ,x,y) :
#     f = interpolation_Cubique
#     return f(f(p00,p01,p02,p03,y),
#              f(p10,p11,p12,p13,y),
#              f(p20,p21,p22,p23,y),
#              f(p30,p31,p32,p33,y),x)



# parameters -------------------------------------------------------------

xwidth  = 10.
xnb_pts = 8

ywidth  = 5.
ynb_pts = 2

X_old = np.linspace(0,xwidth,xnb_pts+1)
Y_old = np.linspace(0,ywidth,ynb_pts+1)

data_nb_pts = (xnb_pts+1)*(ynb_pts+1)
data = rd.normal(1,0.1,data_nb_pts).reshape(ynb_pts+1,xnb_pts+1)

X_new = np.linspace(0,xwidth,xnb_pts*2+1)
Y_new = Y_old
#Y_new = np.linspace(0,ywidth,ynb_pts*2+1)

#print X_old
#print X_new

# script -----------------------------------------------------------------

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")

for j in range(len(Y_new)) : 
    interpol_data = np.zeros(data.shape)
    interpol_1d_x = np.zeros(len(X_new))
    A,B,C,D = 0.,0.,0.,0.
    i_old = 0
    #j_old = ??
    # uncomment this when ready to dev border conditions
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #X0,X1,X2,X3 = x   [i_old-2:i_old+2]
    #P0,P1,P2,P3 = data[i_old-2:i_old+2]
    #A,B,C,D     = update_coefficients(X0,P0,X1,P1,X2,P2,X3,P3)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for i in range(len(X_new)) :
        update_required = False
        while (X_old[i_old] < X_new[i]) :
            i_old += 1
            update_required = True
        # print "j,i,i_old =", j,i,i_old
        # print "\t X_old[i_old] = ",X_old[i_old], "X_new[i] =", X_new[i]
        if update_required and i_old > 1 and i_old < xnb_pts :
            X0,X1,X2,X3 = X_old [i_old-2:i_old+2]
            P0,P1,P2,P3 = data  [j,i_old-2:i_old+2]
            A,B,C,D = update_coefficients(X0,P0,X1,P1,X2,P2,X3,P3)
        if i_old > 1 and i_old < xnb_pts*2 : #forget about the borders for now
            interpol_1d_x[i] = interpolation_Cubique(A,B,C,D,X_new[i])

        #print interpol_1d_x
    ax.scatter(X_old,Y_old[j]*np.ones(len(X_old)),data[j])
    ax.plot(X_new,Y_new[j]*np.ones(len(X_new)),interpol_1d_x, color='r')

plt.ion();plt.show();plt.ioff();raw_input()# uncomment for tests purposes
#plt.savefig("coucou.png")

    
