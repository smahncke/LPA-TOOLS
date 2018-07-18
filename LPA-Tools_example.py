
### IMPORTS CLASSES
import numpy as np
import math
import os
import matplotlib.pyplot as plt
from pylab import title, show
from scipy.constants import c,e
from laser import laser_profiles as las
from elementary_processes import ionization as ion
from elementary_processes import potential as pot
from elementary_processes import program_tools as tool
from particles import particle_data as ptcl
from particles import plasma as plsm


#----------------
# LASER PARAMETER
#----------------

#Peak a0 of  the laser pulse
a0 = 2.1

#Wavelength of the laser (in microns)
lambda_0 = 800e-3

#Waist of the laser in the focus (in microns)
w0 = 22

#Laser duration (in microns)
ctau = 10.1

#Laser focus (in microns)
zf = 50

#-------------------
# PARTICLE PARAMETER
#-------------------

#Max gas density (in 1/cm^3)
ne = 1e24

#The used element (use the shortcut)
element = 'N'

#The energy level that should get ionized
ion_level = 6

#Energy of the ionized electrons after ionization (in Joule)
energy = 0*e


#---------------------
# SIMULATION PARAMETER
#---------------------

#Array over the laser pulse
zz = np.linspace(3*ctau,-3*ctau,1e6)


###################################################


print("")
print('{:*^75}'.format(' Laser plasma accelerator tools '))
print("")

#------------------------------------
# Calculate the ionization parameters
#------------------------------------

#Ionization energy
U_i = ptcl.get_ion_energy(element,ion_level)

#Ionization probability
prob = ion.ionization_probability(zz,a0,w0,ctau,zf,lambda_0,U_i,energy,)

#Ionization degree
degree = ion.ionization_degree(zz,prob)

#Laser field
laser_field = las.gaussian_field(zz, a0,w0,ctau,zf,lambda_0)

#Laser envelope
laser_envelope = las.gaussian_a0(zz, a0,ctau)

#Wakefield
wake = pot.wakefield(zz,a0,ctau,ne)

#Ionization energy distribution
#distribution = ion.ionization_energy_distribution(zz,a0,w0,ctau,zf,lambda_0,U_i,energy_range=(0,5,0.1))

#Critical plasma density
n_c = plsm.get_critical_density(lambda_0*1e-6)

#-------------------
# Print the results
#-------------------

print("")
print('{:-^75}'.format(' Results '))
print("")
print("The "+str(ion_level)+". niveau of the element '"+str(element)+"' has an ionization energy of "+str(U_i)+" eV.")
print("")
print("The critical plasma density for a laser with a wavelength of "+str(lambda_0*1e-6)+" nm")
print("is "+str(n_c)+" 1/m**3")
print("")
print("With a peak a0 of "+str(a0)+", the max. ionization probability is "+str(100*prob.max())+" %.")
print("")
print("The final degree of ionization behind the laser pulse is "+str(100*degree[len(degree)-1])+" %")
print("if the electrons should have a final kinetic energy of "+str(energy/e)+" eV")
print("")
print('{:-^75}'.format(' Website '))
print("* Visit https://github.com/smahncke/LPA-TOOLS for documentation *")
print("")

#-----------------
# Plot the results
#-----------------

tool.plotter(zz,[degree,laser_field,wake],x_label="z",y_label="[a.u.]",\
				plot_title = r"LPA-Tools example: $a_0 = $"+str(a0)+", "+str(element)+str(ion_level)+r"+ $\rightarrow$ "+str(element)+str(ion_level+1)+"+")

#tool.plotter(distribution[0],distribution[1], x_label="Energy [eV]", y_label="Degree [%]", plot_title="Ionization energy distribution")



