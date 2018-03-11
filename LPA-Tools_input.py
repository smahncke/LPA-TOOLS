
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

a0 = 1.5
lambda_0 = 800e-3
w0 = 22
ctau = 10.1
zf = 50
r = 0
t = 0

#-------------------
# PARTICLE PARAMETER
#-------------------

element = 'N'
ion_level = 7
energy = 0


#---------------------
# SIMULATION PARAMETER
#---------------------

zz = np.linspace(3*ctau,-3*ctau,1e6)



#--------
# RESULTS
#--------

#Ionization energy
U_i = ptcl.get_ion_energy(element,ion_level)

#Ionization probability
prob = ion.ionization_probability(U_i,energy,a0,zz,r,t,w0,ctau,zf,lambda_0)

#Ionization degree
degree = ion.ionization_degree(prob,zz)

#Laser field
laser_field = las.gaussian_field(a0,zz,r,t,w0,ctau,zf,lambda_0)

#Laser envelope
laser_envelope = las.gaussian_envelope(a0,zz,ctau)






