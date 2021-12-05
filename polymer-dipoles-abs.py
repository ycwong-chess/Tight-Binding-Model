import numpy as np                              # include the built-in numerical methods library for Python
import matplotlib.pyplot as pl                  # include the built-in plotting library
from functions import * 
from conformer import *
                
colors = ['green', 'yellow', 'orange', 'red', 'purple']

#Hello

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
loud = 1                              # flag to determine how dipoles are chosen
N_o = 4                               # number of oligomers

oligomer_list = [Conformer(Nx, Alpha, i+1) for i in range(N_o)]
EList = np.zeros([N_o, Nx])
tdmList = np.zeros([N_o, Nx, 3])
abs_set = np.zeros([N_o, NE])         # set of absorption arrays
lum_set = np.zeros([N_o, NE])  
    
'''Set up site positions and dipoles'''
for i in range(N_o):
    oligomer_list[i].formation()
    H = Ham(oligomer_list[i].DipoleList)
    print('Coupled Hamiltonian is set up:\n')                  
    evals, evecs = np.linalg.eigh(H)
    
    for j in range(0, Nx):
        if (loud):
#            print('State energy', evals[j], evecs[j])
            tdm = (0,0,0)
        for k in range (0, Nx):
            tdm = tdm + oligomer_list[i].DipoleList[k].dir * evecs[k, j]
        tdmList[i][j][:] = tdm
        EList[i][j] = evals[j]
        
#    for j in range (0, Nx):
#        tdm = tdmList[i][j]
#        tdm2 = np.dot(tdmList[i][j],tdmList[i][j])
#        if (loud):
#            print('State ',j,'; Energy', EList[i][j], '; tdm', tdm,'; magnitude',"{:6.4f}".format(tdm2))
        
    '''Set up arrays for plotting'''
    # plot absorption coefficient spectrum
    Eplot, absplot = Plot_Absorption(tdmList[i], EList[i], Emin, Emax, NE)
    abs_set[i] = absplot
    lumplot = Luminescence(tdmList[i], EList[i], Emin, Emax, NE)
    lum_set[i] = lumplot
    

'''Plot average absorption'''
for i in range(N_o):
    pl.plot(Eplot, abs_set[i], color = colors[i], alpha = 0.5)
average_absorption = np.average(abs_set, axis = 0)
pl.plot(Eplot, average_absorption, color='b', lw = 2)
pl.ylim(0, max(average_absorption))
pl.xlabel("Photon energy / eV")
pl.ylabel("Average absorption Coefficient")
pl.title('Average absorption spectrum for %d oligomers of length %d' % (N_o, Nx))
#print("Modelled absorption for ", Nx,'sites and',Ntrap,' traps at angle',"{:6.2f}".format(Alpha*180/np.pi),'degrees')
pl.show()
    
'''Plot average luminescence'''
average_luminescence = np.average(lum_set, axis = 0)
pl.plot(Eplot, average_luminescence, color='r', lw=2)
pl.ylim(0, max(average_luminescence))
pl.yscale('log')
pl.ylim(0.1, 10000)
pl.xlabel("Photon energy / eV")
pl.ylabel("Luminescence")
pl.title('Modelled luminescence for %d oligomers of length %d' % (N_o, Nx))
#print("Modelled luminescence for ", Nx,'sites and', Ntrap,' traps at angle',"{:6.2f}".format(Alpha*180/np.pi),'degrees')
pl.show()




