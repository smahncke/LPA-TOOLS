
### IMPORTS CLASSES
import numpy as np
import math
import os
import matplotlib.pyplot as plt
from pylab import title, show
from scipy.constants import c
from laser import laser_profiles as las
from elementary_processes import ionization as ion
from particles import particle_data as ptcl


#----------------
# LASER PARAMETER
#----------------

#Peak a0 of  the laser pulse
a0 = 1.5

#Wavelength of the laser (in microns)
lambda_0 = 800e-3

#Waist of the laser in the focus (in microns)
w0 = 22

#Laser duration (in femtoseconds)
ctau = 10.1

#Laser focus (in microns)
zf = 50

#-------------------
# PARTICLE PARAMETER
#-------------------

#The used element (use the shortcut)
element = 'N'

#The energy level that should get ionized
ion_level = 7

#Energy of the ionized electrons after ionization (in MeV)
energy = 0


#---------------------
# SIMULATION PARAMETER
#---------------------

#Array over the laser pulse
zz = np.linspace(3*ctau,-3*ctau,1e6)

###################################################


#------------------------------------
# Calculate the ionization parameters
#------------------------------------

#Ionization energy
U_i = ptcl.get_ion_energy(element,ion_level)

#Ionization probability
prob = ion.ionization_probability(U_i,energy,a0,zz,w0,ctau,zf,lambda_0)

#Ionization degree
degree = ion.ionization_degree(prob,zz)

#Laser field
laser_field = las.gaussian_field(a0,zz,w0,ctau,zf,lambda_0)

#Laser envelope
laser_envelope = las.gaussian_envelope(a0,zz,ctau)


#-------------------
# Prints the results
#-------------------

print("----------------------------------------------------------")
print("Results of the LPA-Tools ionization calculation script:")
print("")
print("The "+str(ion_level)+". niveau of the element '"+str(element)+"' has an ionization energy of "+str(U_i)+" MeV.")
print("")
print("With a peak a0 of "+str(a0)+", the max. ionization probability is "+str(100*prob.max())+" %.")
print("")
print("The final degree of ionization behind the laser pulse is "+str(100*degree[len(degree)-1])+" %.")
print("----------------------------------------------------------")






