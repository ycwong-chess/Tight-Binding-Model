#!/usr/bin/env python

import numpy as np                              # include the built-in numerical methods library for Python
import matplotlib.pyplot as pl                  # include the built-in plotting library
import random as ran # import built-in random number library
from scipy.spatial.transform import Rotation as Rot

class Dipole:
    '''Theta and phi are spherical coordinates'''
    def __init__(self, x=0, y=0, z=0, E=0, theta=0, phi=0, dstrength=1):
        self.x=x
        self.y=y
        self.z=z
        self.E=E0
        self.theta=theta
        self.phi=phi
        self.dstrength=dstrength
        self.dir=np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])*dstrength
#    @classmethod
    def setE(cls, Eval):
        cls.E=Eval
#        print('in setE!', Eval)


def Sep(D1,D2):
    return np.sqrt((D1.x-D2.x)**2+(D1.y-D2.y)**2+(D1.z-D2.z)**2)

def SepVec(D1,D2):
    SepV=np.array([D1.x-D2.x, D1.y-D2.y, D1.z-D2.z])*(1/Sep(D1,D2))
    return SepV

'''Can be optimised''' 
'''No superexchange component'''
def Ham(DList):
    N=len(DList)
    H=np.zeros((N,N))
    for i in range (0,N):
        H[i,i]=DList[i].E
    for i in range (0,N):
        for j in range (0,i):
            R=Sep(DList[i],DList[j])
            RVec=SepVec(DList[i],DList[j])
            H[i,j]=np.dot(DList[i].dir,DList[j].dir)/(R**3)-(3*np.dot(DList[i].dir,RVec)*np.dot(DList[j].dir,RVec))/(R**3)
            H[j,i]=H[i,j]
    return H

'''Re-write'''
'''That's just a unit matrix'''
def Ham0(DList):
    N=len(DList)
    H=np.zeros((N,N))
    for i in range (0,N):
        H[i,i]=DList[i].E
    for i in range (0,N):
        for j in range (0,i):
            H[i,j]=0
            H[j,i]=0
    return H

def Gauss(x, x0,sig):
    return 1/(sig*np.sqrt(2*np.pi))*np.exp(-((x-x0)**2/(2*sig*sig)))

def Plot_Dipoles(DList):
#	"Function to plot the dipoles"
    xmin=0.
    xmax=0.
    ymin=0.
    ymax=0.
    N=len(DList)
    for i in range(0,N):
        dx=(DList[i].dstrength)*np.sin(DList[i].theta)*np.cos(DList[i].phi)
        dy=(DList[i].dstrength)*np.sin(DList[i].theta)*np.sin(DList[i].phi)
#        print('dx dy',dx,dy,DList[i].dstrength)
        if DList[i].E<0.95*E0:
            col='r'
        else:
            col='g'
        pl.plot(DList[i].x,DList[i].y,color=col,marker="o", alpha=0.5, markersize=20,ls='None')
        pl.arrow(DList[i].x-dx/2,DList[i].y-dy/2,dx,dy,width=.01)
        pl.text(DList[i].x,DList[i].y, '%d' % i)
#        print('arrow',i,'start finish',DList[i].x-dx/2,DList[i].y-dy/2,DList[i].x+dx/2,DList[i].y+dy/2)
        '''Do not need the 2 lines below?'''
        dx=DList[i].dstrength*np.sin(DList[i].theta)*np.cos(DList[i].phi)
        dy=DList[i].dstrength*np.sin(DList[i].theta)*np.cos(DList[i].phi)
        xmin=np.minimum(DList[i].x,xmin)
        xmax=np.maximum(DList[i].x,xmax)
        ymin=np.minimum(DList[i].y,ymin)
        ymax=np.maximum(DList[i].y,ymax)
#        print('limits',xmin, xmax, ymin, ymax)
        
#        pl.arrow()
    pl.xlim(xmin,xmax)
    pl.ylim(ymin,ymax)
#    pl.Axes.axes.set_aspect()
#	pl.title("dipoles"%(MO_num))
    pl.show()

def Plot_Absorption(tList,min,max,NE):
    dE=(max-min)/NE
    Eplot = []
    absplot = []
    N=len(tList)
    for iE in range (0, NE):
        E=Emin+iE*dE
        absn=0
        for s in range (0,N):
            absn+=E*np.dot(tdmList[s],tdmList[s])*Gauss(E,EList[s],Esig)
        Eplot.append(E)
        absplot.append(absn/N)
    pl.plot(Eplot,absplot,color='b',lw=2, label='Coupled system')
    pl.ylim(0,np.amax(absplot))
    pl.xlabel("Photon energy / eV")
    pl.ylabel("Absorption Coefficient")
    pl.legend()
    print("Modelled absorption for ", Nx,'sites and',Ntrap,' traps at angle',"{:6.2f}".format(Alpha*180/np.pi),'degrees')
    pl.show()
  
#    fig=pl.figure()
#    ax=fig.add_subplot(111,projection='3d')			# Initialise 3D plot

# main program
# Define parameters
E0=1.8                                # energy of exciton before disorder
D0=.1                                # coefficient of diple coupling strenth
kT=0.025                            # kT
Esig=0.01                           # broadening of transitions in spectra
Nx=20                            # number of monomers in chain
Alpha=np.pi/12                      # angle between long axes of consecutive monomers
Length=1                            # Length of a monomer
Ntrap=0
Emin=E0-0.1
Emax=E0+0.1
NE=100
#N=4                                 # number of dipoles in the system
loud=False                              # flag to determine how dipoles are chosen

#D1 = Dipole(1,0,0,E0,0,0,D0)
#print(D1.x, D1.y, D1.z)
#print(D1.dir)
#print(D1.E)

DipoleList = []
EdgeList = []
zaxis=np.array([0,0,1])
#for i in range (0,Nx):
#    DipoleList.append(Dipole(i,j,0,E0,np.arccos(2*ran.random()-1),2*np.pi*ran.random(),D0))
#        DipoleList.append(Dipole(i,j,0,E0,np.pi/2,2*np.pi*ran.random(),D0))
#       DipoleList.append(Dipole(i,j,0,E0,np.pi/2,0,D0))
N=len(DipoleList)
print('length of list =',N)
#v=DipoleList[0].dir

#initialise list of dipoles

v0=np.array([0,0,0])
v=np.array([1,0,0])
n=([0,0,1])
k=np.cross(v,n)
print('Vector before rotation', v,'n',n,'k',k)
#print (v[0],v[1],v[2])
vtheta=np.arccos(v[2])
vphi=np.arctan2(v[1],v[0])
EndPoint=v0+Length*v
MidPoint=v0+Length*v*0.5
#EdgeList.append(EndPoint)
EdgeList.append(EndPoint)
DipoleList.append(Dipole(MidPoint[0],MidPoint[1],MidPoint[2],E0,vtheta,vphi,D0))
# enter loop having set up vectors v, n, k
for i in range (0,Nx):
#    print('v.n',np.dot(v,n),'vxn',np.cross(v,n))
    r=Rot.from_rotvec(Alpha*k)
    v2=r.apply(v)
    n=r.apply(n)
#    print('Vector after lift', v2)
    Angle=ran.random()*np.pi*2
    RotVec=np.transpose(v)
    r=Rot.from_rotvec(Angle*RotVec)
    v=r.apply(v2)
    n=r.apply(n)
# now rotate the normal by a random angle in 0 to 2pi
    Angle=ran.random()*np.pi*2
    r=Rot.from_rotvec(Angle*v)
    n=r.apply(n)
#    print('v.n',np.dot(v,n))
    k=np.cross(v,n)
#    print('Vector after rotation', v,'|v|',np.linalg.norm(v),'n',n,'k',k,'|k|',np.linalg.norm(k))
    vtheta=np.arccos(v[2])
    vphi=np.arctan2(v[1],v[0])
    EndPoint=EdgeList[-1]+Length*v
    MidPoint=EdgeList[-1]+Length*v*0.5
    EdgeList.append(EndPoint)
    DipoleList.append(Dipole(MidPoint[0],MidPoint[1],MidPoint[2],E0,vtheta,vphi,D0))
#    print('point',EdgeList[-1],'vec',v,'length',np.linalg.norm(v))    
#    print('point',EndPoint,'v',v)
print('length of points list',len(EdgeList),'length of dipole list',len(DipoleList))
for i in range (0,Nx):
    v=EdgeList[i+1]-EdgeList[i]
    print(i,'point',EdgeList[i],'vec',v,'length',np.linalg.norm(v),'dipole', DipoleList[i].dir)    
#print(i,'point',EdgeList[-1],'vec',v,'length',np.linalg.norm(v),'dipole', DipoleList[-1].dir)    
#RotVec=k
#print('Rotation Vector', RotVec, 'Angle',Angle)
#r=Rot.from_rotvec(Angle*RotVec)
#v2=r.apply(v)
#print('Vector before rotation', v)
#print('Vector after lift', v2)
#Angle=ran.random()*np.pi*2
#RotVec=np.transpose(v)
#print('v',v,'Rotation Vector', RotVec, 'Angle',Angle)
#r=Rot.from_rotvec(Angle*RotVec)
#v2=r.apply(v2)
#print('Vector after rotation', v2)

print('dipole positions set')
N=len(DipoleList)
for i in range (0,len(DipoleList)):
    dip=DipoleList[i]
    print('no.',i,' pos. ',dip.x,dip.y,dip.z, 'dipole:',dip.dir)
Plot_Dipoles(DipoleList)

for i in range (0,N):
    for j in range (0,i):
        Dip1=DipoleList[i]
        Dip2=DipoleList[j]
        R=Sep(DipoleList[i],DipoleList[j])
#        print(i,Dip1.x,Dip1.y,Dip1.z,'|',j,Dip2.x,Dip2.y,Dip2.z,'|',SepVec(Dip1,Dip2))
#        print('separation of elements',i,' and ',j,' is ',"{:6.3f}".format(R))



# set up site positions and dipoles (and energies?)
H=Ham(DipoleList)
#print('Ham set up:\n',H)            
print('Coupled Ham set up:\n')            

tdmList = []
EList =[]
tdm=np.array([0, 0, 0])
evals,evecs=np.linalg.eigh(H)
for s in range(0,N):
    if (loud):
        print('State energy',evals[s], evecs[s])
    tdm=(0,0,0)
    for i in range (0,N):
        tdm=tdm+DipoleList[i].dir*evecs[i,s]
#    print('State tdm:',tdm)
    tdmList.append(tdm)
    EList.append(evals[s])
    
for s in range (0,N):
    tdm=tdmList[s]
    tdm2=np.dot(tdmList[s],tdmList[s])
    if (loud):
        print('State ',s,'Energy', EList[s],tdm,'mag',"{:6.4f}".format(tdm2))
    

# set up graphs to plot

# plot absorption coefficient spectrum
Plot_Absorption(tdmList,Emin,Emax,NE)
# Set up and plot luminescnce
dE=(Emax-Emin)/NE
Eplot = []
lumplot = []
for iE in range (0, NE):
    E=Emin+iE*dE
    lum=0
    for s in range (0,N):
        lum+=E**3*np.exp(-(EList[s]-E0)/kT)*np.dot(tdmList[s],tdmList[s])*Gauss(E,EList[s],Esig)
    Eplot.append(E)
    lumplot.append(lum)
pl.plot(Eplot,lumplot,color='r',lw=2, label='Coupled system')
pl.ylim(0,np.amax(lumplot))
pl.yscale('log')
pl.ylim(0.1,10000)
pl.xlabel("Photon energy / eV")
pl.ylabel("Luminescence")
pl.legend()
#pl.title("Modelled luminescence for ",%Nx,'sites with alpha',%Alpha)
print("Modelled luminescence for ", Nx,'sites and',Ntrap,' traps at angle',"{:6.2f}".format(Alpha*180/np.pi),'degrees')
pl.show()

#Ham=ring_H(N,E0,t)					# Call function to construct an N-length chain

#evals,evecs=np.linalg.eigh(Ham)		# Solve H\psi = E\psi

#xs=np.arange(0,N)					# Get x and y positions for orbitals (y set to zero)
#ys=np.zeros(N)

#xc,yc= Polygon(N,1,0,0)
#for MO_num in range(0,N):
   
#    Plot_orbitals(evecs,N,xc,yc,MO_num)	# Plot orbitals of the chains
#    print ("Chain eigenvalue: ",MO_num, evals[MO_num])


