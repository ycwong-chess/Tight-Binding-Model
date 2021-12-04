#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 11:56:19 2021

@author: brilliant
"""

import numpy as np
import matplotlib.pyplot as pl                  # include the built-in plotting library
import random as ran                            # import built-in random number library
from scipy.spatial.transform import Rotation as Rot
from functions import *

E0 = 1.8                              # energy of exciton before disorder
D0 = 0.1                              # coefficient of diple coupling strenth
Length = 1                            # Length of a monomer

class Dipole:
    def __init__(self, x=0, y=0, z=0, E=0, theta=0, phi=0, dstrength=1):
        self.x = x
        self.y = y
        self.z = z
        self.E = E0
        self.theta = theta #spherical coord
        self.phi = phi
        self.dstrength = dstrength
        self.dir = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])*dstrength
        
    def setE(self, Eval):
        self.E = Eval
        
class Conformer: 
    def __init__(self, number, angle, index):
        self.units = number
        self.alpha = angle
        self.index = index
        self.DipoleList = []
        self.EdgeList = []
#         Set up vectors v, n, k
        v0 = np.array([0,0,0])
        self.v = np.array([1,0,0])
        self.n = np.array([0,0,1])
        self.k = np.cross(self.v, self.n)
        
        self.vtheta = np.arccos(self.v[2])
        self.vphi=np.arctan2(self.v[1],self.v[0])
        
        self.EndPoint = v0 + Length * self.v
        self.MidPoint = v0 + Length * self.v/2
        self.EdgeList.append(self.EndPoint)
        self.DipoleList.append(Dipole(self.MidPoint[0], self.MidPoint[1], self.MidPoint[2], \
                                      E0, self.vtheta, self.vphi, D0))
        
    def formation(self):
        '''Form a list of dipoles'''       
        for i in range (0, self.units):
#           Rotation about k by alpha
            r = Rot.from_rotvec(self.alpha * self.k)
            v2 = r.apply(self.v)
            self.n = r.apply(self.n)
#           Rotation about v by random between 0 & 2pi
            Angle1 = ran.random()*np.pi*2
            r = Rot.from_rotvec(Angle1 * self.v) #*RotVec instead of v; RotVec = np.transpose(v)
            self.v = r.apply(v2)
            self.n = r.apply(self.n)
#           Rotating the normal about the new v by random between 0 & 2pi
            Angle2 = ran.random()*np.pi*2
            r = Rot.from_rotvec(Angle2 * self.v)
            self.n = r.apply(self.n)
#           Updating the rest of the variables
            self.k = np.cross(self.v,self.n)
            self.vtheta = np.arccos(self.v[2])
            self.vphi = np.arctan2(self.v[1], self.v[0])
            
            self.EndPoint = self.EdgeList[-1] + Length * self.v
            self.MidPoint = self.EdgeList[-1] + Length * self.v/2
            self.EdgeList.append(self.EndPoint)
            self.DipoleList.append(Dipole(self.MidPoint[0], self.MidPoint[1], self.MidPoint[2],\
                                       E0, self.vtheta, self.vphi, D0))
#        print('Dipole positions set (length of list: %d)' % len(self.DipoleList))
        Plot_Dipoles(self.DipoleList, self.index)



