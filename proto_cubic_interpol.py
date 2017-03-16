#!/usr/bin/python
#-*-coding:utf-8-*-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# written march 2017
# Clément Robert
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# imports ----------------------------------------------------------------

import numpy as np
from numpy import random as rd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# defintions -------------------------------------------------------------

def update_coefficients(x0,p0,x1,p1,x2,p2,x3,p3,a,b,c,d) :
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
    return a,b,c,d # this is the pythonic way to give back the values. Not useful in C


def third_degree_polynom(a,b,c,d,x) :
    return a*x**3 + b*x**2 + c*x + d

tdp = third_degree_polynom #alias


def update_indexes(i_old,j_old,xmax,ymax) :#todo : fake C incorporate returns as arguments


    #dev note : those booleans should be global variables
    useXghostNEG1=useXghost0=useXghost1=False



    # handling exceptions near azimuthal "borders" (periodic, really)
    if j_old == 1 :                    # case 1
        seed        = 0
        s0,s1,s2,s3 = xnb_pts-1,seed,seed+1,seed+2
        useXghostNEG1 = True
    elif j_old == xnb_pts-1 :          # case 2
        seed = xnb_pts-3
        s0,s1,s2,s3 = seed,seed+1,seed+2,0
        useXghost0 = True
    elif X_new[j] > X_old[xnb_pts-1] : # case 3
        #dev note : this condition should be better written using % [2PI]
        seed = xnb_pts-2
        s0,s1,s2,s3 = seed,seed+1,0,1
        useXghost0=useXghost1 = True
    else :                             # general case
        seed = j_old-2
        s0,s1,s2,s3 = seed,seed+1,seed+2,seed+3


    # handling exceptions near radial borders
    if i_old == 0 :           # case 1
        shift = 2
    elif i_old == 1 :         # case 2
        shift = 1
    elif i_old == ynb_pts-1 : # case 3
        shift = -1
    elif i_old == ynb_pts :   # case 4
        shift = -2
    else :                    # general case
        shift = 0

    i0 = i_old-2+shift
    # computation at long last
    l         = s0 + i0*xnb_pts
    ljp       = s1 + i0*xnb_pts
    ljpp      = s2 + i0*xnb_pts
    ljppp     = s3 + i0*xnb_pts

    lip       = l     +   xnb_pts
    lipp      = l     + 2*xnb_pts
    lippp     = l     + 3*xnb_pts

    lipjp     = ljp   +   xnb_pts
    lippjp    = ljp   + 2*xnb_pts
    lipppjp   = ljp   + 3*xnb_pts

    lipjpp    = ljpp  +   xnb_pts
    lippjpp   = ljpp  + 2*xnb_pts
    lipppjpp  = ljpp  + 3*xnb_pts

    lipjppp   = ljppp +   xnb_pts
    lippjppp  = ljppp + 2*xnb_pts
    lipppjppp = ljppp + 3*xnb_pts

    return s0,s1,s2,s3,l,\
        lip,    lipp,    lippp,\
        ljp,    ljpp,    ljppp,\
        lipjp,  lipjpp,  lipjppp,\
        lippjp, lippjpp, lippjppp,\
        lipppjp,lipppjpp,lipppjppp,\
        useXghostNEG1,useXghost0,useXghost1,\
        i0,shift


# parameters -------------------------------------------------------------

xwidth  = 2.*np.pi
xnb_pts = 7
xwidth_old = xwidth*(1.-1./xnb_pts)

ywidth  = 10.
ynb_pts = 8
ystep_old = ywidth / (ynb_pts +2)

X_old = np.linspace(0.,xwidth_old,xnb_pts)
Y_old = np.linspace(ystep_old,ywidth-ystep_old,ynb_pts)

XMAX = X_old[xnb_pts-2]
YMAX = Y_old[ynb_pts-2]


#assuming that xmin=0, and xwidth=2Pi, which should always be the case
X_GHOSTneg1 = X_old[-1] - 2*np.pi
X_GHOST0    = X_old[0 ] + 2*np.pi
X_GHOST1    = X_old[1 ] + 2*np.pi

data_nb_pts = xnb_pts*ynb_pts
data1d = rd.normal(1,0.1,data_nb_pts)

x_enhance_factor = 10
y_enhance_factor = 10
xwidth_new = xwidth*(1.-1./(xnb_pts*x_enhance_factor-1))
X_new = np.linspace(0,xwidth_new,xnb_pts*x_enhance_factor-4)
Y_new = np.linspace(0,ywidth,ynb_pts*y_enhance_factor-1)


# script -----------------------------------------------------------------

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")

Y00=Y01=Y02=Y03= \
Y10=Y11=Y12=Y13= \
Y20=Y21=Y22=Y23= \
Y30=Y31=Y32=Y33= 0.0


AA=BB=CC=DD=0.
A0=B0=C0=D0=0.
A1=B1=C1=D1=0.
A2=B2=C2=D2=0.
A3=B3=C3=D3=0.
A =B =C =D =0.


for i in range(len(Y_new)) :
    y_new = Y_new[i]
    # (re)init interpolated data
    interpol_data = np.zeros(len(X_new))
    interpol_1d_x = np.zeros(len(X_new))

    # update y index
    i_old = 0
    while (Y_old[i_old] < Y_new[i]) :
        i_old += 1
        if i_old == ynb_pts :
            break


    for j in range(len(X_new)) :
        x_new = X_new[j]
        y_new_dummy = y_new
        update_required = False

        # update x index
        j_old = 0
        while (X_old[j_old] <= X_new[j]) :
            j_old += 1
            update_required = True
            if j_old == xnb_pts :
                break

        s0,s1,s2,s3,l,\
        lip,    lipp,    lippp,\
        ljp,    ljpp,    ljppp,\
        lipjp,  lipjpp,  lipjppp,\
        lippjp, lippjpp, lippjppp,\
        lipppjp,lipppjpp,lipppjppp,\
        useXghostNEG1,useXghost0,useXghost1,\
        i0,shift\
        = update_indexes(i_old,j_old,XMAX,YMAX)

        if update_required :
            # the routine itself might be written in all generality but in practice
            # we know we will only use evenly spaced grids IN THETA (aka x here)
            X00,X01,X02,X03 = \
            X10,X11,X12,X13 = \
            X20,X21,X22,X23 = \
            X30,X31,X32,X33 = X_old [s0], X_old [s1], X_old [s2], X_old [s3]

            if   useXghostNEG1                 : #case 1
                X00=X10=X20=X30 = X_GHOSTneg1
            elif useXghost0 and not useXghost1 : #case 2
                X03=X13=X23=X33 = X_GHOST0
            elif useXghost0 and     useXghost1 : #case 3
                X02=X12=X22=X32 = X_GHOST0
                X03=X13=X23=X33 = X_GHOST1

            # however, we don't enconter any concerning issue with log-spaced
            # radial grids but the syntax ought to be different
            Y00=Y01=Y02=Y03 = Y_old [i0  ]
            Y10=Y11=Y12=Y13 = Y_old [i0+1]
            Y20=Y21=Y22=Y23 = Y_old [i0+2]
            Y30=Y31=Y32=Y33 = Y_old [i0+3]

            if shift == 2 :
                y_new_dummy = Y00
            elif shift == -2 :
                y_new_dummy = Y30

            # field values
            P00,P01,P02,P03 = data1d[l],    data1d[ljp],    data1d[ljpp],     data1d[ljppp]
            P10,P11,P12,P13 = data1d[lip],  data1d[lipjp],  data1d[lipjpp],   data1d[lipjppp]
            P20,P21,P22,P23 = data1d[lipp], data1d[lippjp], data1d[lippjpp],  data1d[lippjppp]
            P30,P31,P32,P33 = data1d[lippp],data1d[lipppjp],data1d[lipppjpp], data1d[lipppjppp]

            # finally, update all the coefficients
            A0,B0,C0,D0 = update_coefficients(X00,P00,X01,P01,X02,P02,X03,P03,A0,B0,C0,D0)
            A1,B1,C1,D1 = update_coefficients(X10,P10,X11,P11,X12,P12,X13,P13,A1,B1,C1,D1)
            A2,B2,C2,D2 = update_coefficients(X20,P20,X21,P21,X22,P22,X23,P23,A2,B2,C2,D2)
            A3,B3,C3,D3 = update_coefficients(X30,P30,X31,P31,X32,P32,X33,P33,A3,B3,C3,D3)

            interpol_1d_x[j] = tdp(AA,BB,CC,DD,x_new)

        # this is where magic happens
        P0 = tdp(A0,B0,C0,D0,x_new)
        P1 = tdp(A1,B1,C1,D1,x_new)
        P2 = tdp(A2,B2,C2,D2,x_new)
        P3 = tdp(A3,B3,C3,D3,x_new)

        A,B,C,D = update_coefficients(Y00,P0,Y10,P1,Y20,P2,Y30,P3,A,B,C,D)
        interpol_value = tdp(A,B,C,D,y_new_dummy)
        interpol_data[j] = interpol_value

        # those lines are used to keep track and debug the process -------------
        #dummy_x = x_new*np.ones(100)
        #dummy_y = np.linspace(0.,ywidth,100)
        #ax.plot(dummy_x,dummy_y,tdp(A,B,C,D,dummy_y),color = 'k', ls = '--')
        # ----------------------------------------------------------------------

    if i_old < ynb_pts :
        ax.scatter(X_old,Y_old[i_old]*np.ones(len(X_old)),data1d[i_old*xnb_pts : (i_old+1)*xnb_pts])
    ax.plot(X_new,y_new*np.ones(len(X_new)),interpol_data, color='c', lw=2, ls='-')

plt.ion();plt.show();plt.ioff();raw_input("press anykey to quit    ")
