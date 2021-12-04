#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 12:38:56 2021

@author: brilliant
"""
import numpy as np                              # include the built-in numerical methods library for Python
import matplotlib.pyplot as pl                  # include the built-in plotting library

E0 = 1.8                              # energy of exciton before disorder
D0 = 0.1                              # coefficient of diple coupling strenth
kT = 0.025                            # kT
Esig = 0.01                           # broadening of transitions in spectra
Nx = 20                               # number of monomers in chain
Alpha = np.pi/12                      # angle between long axes of consecutive monomers
Length = 1                            # Length of a monomer
Ntrap = 0
Emin = E0-0.1
Emax = E0+0.1
NE = 100
loud = 0                              # flag to determine how dipoles are chosen
N_o = 4 

def Sep(D1,D2):
    return np.sqrt((D1.x-D2.x)**2+(D1.y-D2.y)**2+(D1.z-D2.z)**2)

def SepVec(D1,D2):
    return np.array([D1.x-D2.x, D1.y-D2.y, D1.z-D2.z])/Sep(D1,D2)

def Ham0(Dlist):
    N = len(Dlist)
    I = np.identity(N)
    return I * np.array([dipole.E for dipole in Dlist])

'''No superexchange component'''
def Ham(DList):
    N=len(DList)
    H = Ham0(DList)
    for i in range (0,N):
        for j in range (0,i):
            R = Sep(DList[i],DList[j])
            RVec = SepVec(DList[i],DList[j])
            H[i,j] = (np.dot(DList[i].dir,DList[j].dir) - \
                    3*np.dot(DList[i].dir,RVec)*np.dot(DList[j].dir,RVec))/(R**3)
            H[j,i] = H[i,j]
    return H

def Gauss(x, x0,sig):
    return 1/(sig*np.sqrt(2*np.pi))*np.exp(-((x-x0)**2/(2*sig**2)))

def Plot_Dipoles(DList, index):
    """Function to plot the dipoles"""
    xmin=0.
    xmax=0.
    ymin=0.
    ymax=0.
    N=len(DList)
    for i in range(0,N):
        dx=(DList[i].dstrength)*np.sin(DList[i].theta)*np.cos(DList[i].phi)
        dy=(DList[i].dstrength)*np.sin(DList[i].theta)*np.sin(DList[i].phi)
        if DList[i].E < 0.95 * E0:
            col='r'
        else:
            col='g'
        pl.plot(DList[i].x,DList[i].y,color=col,marker="o", alpha=0.5, markersize=20,ls='None')
        pl.arrow(DList[i].x-dx/2,DList[i].y-dy/2,dx,dy,width=.01)
        pl.text(DList[i].x,DList[i].y, '%d' % i)
        pl.title('Conformer %d' % (index))
        pl.xlabel('X coordinate')
        pl.ylabel('Y coordinate')
        xmin=np.minimum(DList[i].x,xmin)
        xmax=np.maximum(DList[i].x,xmax)
        ymin=np.minimum(DList[i].y,ymin)
        ymax=np.maximum(DList[i].y,ymax)
    pl.xlim(xmin,xmax)
    pl.ylim(ymin,ymax)
#   pl.Axes.axes.set_aspect()
#	pl.title("dipoles"%(MO_num))
    pl.show()

def Plot_Absorption(tdmList, Elist,min,max,NE):
    dE=(max-min)/NE
    Eplot = []
    absplot = []
    N = len(tdmList)
    for iE in range (0, NE):
        E = min + iE*dE
        absn=0
        for s in range (0, N):
            absn += E * np.dot(tdmList[s],tdmList[s]) * Gauss(E,Elist[s], Esig)
        Eplot.append(E)
        absplot.append(absn/N)
#    pl.plot(Eplot,absplot,color='b',lw=2, label='Coupled system')
#    pl.ylim(0,np.amax(absplot))
#    pl.xlabel("Photon energy / eV")
#    pl.ylabel("Absorption Coefficient")
#    pl.legend()
##    print("Modelled absorption for ", Nx,'sites and',Ntrap,' traps at angle',"{:6.2f}".format(Alpha*180/np.pi),'degrees')
#    pl.show()
    return Eplot, absplot

def Luminescence(tdmList, Elist, min, max, NE):
    dE=(max-min)/NE
    lumplot = []
    N = len(tdmList)
    for iE in range (0, NE):
        E = min + iE*dE
        lum = 0
        for s in range (0, N):
            lum += E**3 * np.exp(-(Elist[s]-E0)/kT) * np.dot(tdmList[s],tdmList[s]) * Gauss(E,Elist[s],Esig)
        lumplot.append(lum)
    return lumplot

    
#    fig=pl.figure()
#    ax=fig.add_subplot(111,projection='3d')			# Initialise 3D plot