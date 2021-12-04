import numpy as np                              # include the built-in numerical methods library for Python
import matplotlib.pyplot as pl                  # include the built-in plotting library
import random as ran                            # import built-in random number library
from scipy.spatial.transform import Rotation as Rot
            
class Dipole:
    '''Theta and phi are spherical coordinates'''
    def __init__(self, x=0, y=0, z=0, E=0, theta=0, phi=0, dstrength=1):
        self.x = x
        self.y = y
        self.z = z
        self.E = E0
        self.theta = theta
        self.phi = phi
        self.dstrength = dstrength
        self.dir = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])*dstrength
#    @classmethod
    def setE(cls, Eval):
        cls.E = Eval

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

def Plot_Dipoles(DList):
    """Function to plot the dipoles"""
    xmin=0.
    xmax=0.
    ymin=0.
    ymax=0.
    N=len(DList)
    for i in range(0,N):
        dx=(DList[i].dstrength)*np.sin(DList[i].theta)*np.cos(DList[i].phi)
        dy=(DList[i].dstrength)*np.sin(DList[i].theta)*np.sin(DList[i].phi)
        if DList[i].E<0.95*E0:
            col='r'
        else:
            col='g'
        pl.plot(DList[i].x,DList[i].y,color=col,marker="o", alpha=0.5, markersize=20,ls='None')
        pl.arrow(DList[i].x-dx/2,DList[i].y-dy/2,dx,dy,width=.01)
        pl.text(DList[i].x,DList[i].y, '%d' % i)
        xmin=np.minimum(DList[i].x,xmin)
        xmax=np.maximum(DList[i].x,xmax)
        ymin=np.minimum(DList[i].y,ymin)
        ymax=np.maximum(DList[i].y,ymax)
    pl.xlim(xmin,xmax)
    pl.ylim(ymin,ymax)
#   pl.Axes.axes.set_aspect()
#	pl.title("dipoles"%(MO_num))
    pl.show()

def Plot_Absorption(tList,min,max,NE):
    dE=(max-min)/NE
    Eplot = []
    absplot = []
    N=len(tList)
    for iE in range (0, NE):
        E = Emin + iE*dE
        absn=0
        for s in range (0,N):
            absn += E*np.dot(tdmList[s],tdmList[s])*Gauss(E,EList[s],Esig)
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
    
"""Main Body"""

'''Define parameters'''
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
loud=1                          # flag to determine how dipoles are chosen

DipoleList = []
EdgeList = []
zaxis=np.array([0,0,1])

'''Set up vectors v, n, k'''
v0 = np.array([0,0,0])
v = np.array([1,0,0])
n = ([0,0,1])
k = np.cross(v,n)
print('Vector before rotation', v, 'n', n, 'k', k)
vtheta = np.arccos(v[2])
vphi=np.arctan2(v[1],v[0])
EndPoint = v0 + Length*v
MidPoint = v0 + Length*v/2
EdgeList.append(EndPoint)
DipoleList.append(Dipole(MidPoint[0], MidPoint[1], MidPoint[2], E0, vtheta, vphi, D0))

'''Form a list of dipoles'''
for i in range (0, Nx):
#   Rotation about k by alpha
    r = Rot.from_rotvec(Alpha*k)
    v2 = r.apply(v)
    n = r.apply(n)
#   Rotation about v by random between 0 & 2pi
    Angle = ran.random()*np.pi*2
    r = Rot.from_rotvec(Angle*v) #*RotVec instead of v; RotVec = np.transpose(v)
    v = r.apply(v2)
    n = r.apply(n)
#   Rotating the normal about the new v by random between 0 & 2pi
    Angle = ran.random()*np.pi*2
    r = Rot.from_rotvec(Angle*v)
    n = r.apply(n)
#   Updating k
    k = np.cross(v,n)
#   print('Vector after rotation', v,'|v|',np.linalg.norm(v),'n',n,'k',k,'|k|',np.linalg.norm(k))
    vtheta = np.arccos(v[2])
    vphi = np.arctan2(v[1],v[0])
    EndPoint = EdgeList[-1] + Length*v
    MidPoint = EdgeList[-1] + Length*v/2
    EdgeList.append(EndPoint)
    DipoleList.append(Dipole(MidPoint[0], MidPoint[1], MidPoint[2], E0, vtheta, vphi,D0))
print('length of points list', len(EdgeList), 'length of dipole list', len(DipoleList))

#for i in range (0, Nx):
#    v = EdgeList[i+1] - EdgeList[i]
#    print(i, 'point', EdgeList[i], 'vec', v, 'length', np.linalg.norm(v), 'dipole', DipoleList[i].dir)    

print('Dipole positions set')

N = len(DipoleList)
for i in range (0,len(DipoleList)):
    dip = DipoleList[i]
    print('no.', i, ' pos. ', dip.x, dip.y, dip.z, 'dipole:', dip.dir)
Plot_Dipoles(DipoleList)

#for i in range (0, N):
#    for j in range (0, i):
#        Dip1 = DipoleList[i]
#        Dip2 = DipoleList[j]
#        R = Sep(DipoleList[i], DipoleList[j])
#        print(i,Dip1.x,Dip1.y,Dip1.z,'|',j,Dip2.x,Dip2.y,Dip2.z,'|',SepVec(Dip1,Dip2))
#        print('separation of elements',i,' and ',j,' is ',"{:6.3f}".format(R))

# set up site positions and dipoles (and energies?)
H = Ham(DipoleList)
print('Coupled Hamiltonian is set up:\n')            

tdmList = []
EList = []
tdm = np.array([0, 0, 0])
evals,evecs = np.linalg.eigh(H)
for s in range(0,N):
    if (loud):
        print('State energy',evals[s], evecs[s])
    tdm = (0,0,0)
    for i in range (0,N):
        tdm = tdm + DipoleList[i].dir*evecs[i,s]
#    print('State tdm:',tdm)
    tdmList.append(tdm)
    EList.append(evals[s])
    
for s in range (0,N):
    tdm = tdmList[s]
    tdm2 = np.dot(tdmList[s],tdmList[s])
    if (loud):
        print('State ',s,'; Energy', EList[s], '; tdm', tdm,'; magnitude',"{:6.4f}".format(tdm2))
    

# set up graphs to plot

# plot absorption coefficient spectrum
Plot_Absorption(tdmList,Emin,Emax,NE)
# Set up and plot luminescnce
dE=(Emax-Emin)/NE
Eplot = []
lumplot = []
for iE in range (0, NE):
    E = Emin + iE*dE
    lum=0
    for s in range (0,N):
        lum += E**3 * np.exp(-(EList[s]-E0)/kT)*np.dot(tdmList[s],tdmList[s])*Gauss(E,EList[s],Esig)
    Eplot.append(E)
    lumplot.append(lum)
    
pl.plot(Eplot,lumplot,color='r',lw=2, label='Coupled system')
pl.ylim(0,np.amax(lumplot))
pl.yscale('log')
pl.ylim(0.1, 10000)
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



