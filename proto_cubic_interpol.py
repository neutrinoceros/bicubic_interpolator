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

def interpolation_BiCubique(x00, y00, p00,
                            x01, y01, p01,
                            x02, y02, p02,
                            x03, y03, p03,
                            x10, y10, p10,
                            x11, y11, p11,
                            x12, y12, p12,
                            x13, y13, p13,
                            x20, y20, p20,
                            x21, y21, p21,
                            x22, y22, p22,
                            x23, y23, p23,
                            x30, y30, p30,
                            x31, y31, p31,
                            x32, y32, p32,
                            x33, y33, p33,x,y) :

    f = interpolation_Cubique#local alias

    #dev note :
    # it would be much more optimal to update coefficients outside of this function
    A0,B0,C0,D0 = update_coefficients(x00,p00,x01,p01,x02,p02,x03,p03)
    A1,B1,C1,D1 = update_coefficients(x10,p10,x11,p11,x12,p12,x13,p13)
    A2,B2,C2,D2 = update_coefficients(x20,p20,x21,p21,x22,p22,x23,p23)
    A3,B3,C3,D3 = update_coefficients(x30,p30,x31,p31,x32,p32,x33,p33)
    P0 = f(A0,B0,C0,D0,x)
    P1 = f(A1,B1,C1,D1,x)
    P2 = f(A2,B2,C2,D2,x)
    P3 = f(A3,B3,C3,D3,x)
    A,B,C,D = update_coefficients(y00,P0,y10,P1,y20,P2,y30,P3)

    return f(A,B,C,D,y)

# parameters -------------------------------------------------------------

xwidth  = 10.
xnb_pts = 6

ywidth  = 5.
ynb_pts = 6

X_old = np.linspace(0,xwidth,xnb_pts+1)
Y_old = np.linspace(0,ywidth,ynb_pts+1)

data_nb_pts = (xnb_pts+1)*(ynb_pts+1)
data = rd.normal(1,0.1,data_nb_pts).reshape(ynb_pts+1,xnb_pts+1)

X_new = np.linspace(0,xwidth,xnb_pts*4+1)
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


    j_old = 0
    i_old = 0
    # uncomment this when ready to dev border conditions
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #X0,X1,X2,X3 = x   [i_old-2:i_old+2]
    #P0,P1,P2,P3 = data[j_old-2:j_old+2,i_old-2:i_old+2]
    #A,B,C,D     = update_coefficients(X0,P0,X1,P1,X2,P2,X3,P3)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for i in range(len(X_new)) :
        Xupdate_required = False
        while (X_old[i_old] < X_new[i]) :
            i_old += 1
            Xupdate_required = True
        if i_old > 1 and X_new[i] < X_old[xnb_pts-1] :
            if Xupdate_required :
                X0,X1,X2,X3 = X_old [i_old-2:i_old+2]
                P0,P1,P2,P3 = data  [j,i_old-2:i_old+2]
                A,B,C,D = update_coefficients(X0,P0,X1,P1,X2,P2,X3,P3)

            interpol_1d_x[i] = interpolation_Cubique(A,B,C,D,X_new[i])
        else : #forget about the borders for now
            pass
    ax.scatter(X_old,Y_old[j]*np.ones(len(X_old)),data[j])
    ax.plot(X_new,Y_new[j]*np.ones(len(X_new)),interpol_1d_x, color='r')

plt.ion();plt.show();plt.ioff();raw_input("press anykey to quit    ")# uncomment for tests purposes
#plt.savefig("coucou.png")
