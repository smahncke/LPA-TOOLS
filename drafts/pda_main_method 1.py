###############################
"""Questions/To do:
1. Soll am Anfang GAS- oder PLASMADICHTE angegeben werden? Falls GASDICHTE, ist n_e_H2_tot dann IMMER 2*n_tot ?
2. Existiert Limit für am_N2, sodass immer 2*n_tot eingehalten wird? (Andernfalls wird n_H2 < 0)
3. Was ist der gewünschte Output? Gas oder Plasmadichte?"""
###############################


### IMPORTS CLASSES
import numpy as np
import math
import time
import os
import sys
import matplotlib.pyplot as plt
from pylab import title, show
from opmd_viewer.addons.pic import LpaDiagnostics as OpenPMDTimeSeries
from scipy.constants import c

# Global settings

ts = OpenPMDTimeSeries('./diags/hdf5/') #The simulation files which should be imported

# Gas settings
n_tot = 1e18*1e6 #Total Plasma Density
am_N2 = 10 #Amount of nitrogen (in percent)

# Density profiles

FWHM = 300e-6 #Full-width at half-max of the gaussian density profile (in m)
ramp_start = 30e-6 #Start point of the resulting plasma upramp (in m)
ramp_length = 500e-6 #The length of the upramp and the downramp (in m)
plateau = 4000e-6 #The length of the plateau (in m)

# Other settings

inj_thres = 6 #The ionization level of N2, at which the electron injection starts

###################################################################################################################################################################


# Gaussian profile
mu = ramp_start + ramp_length #mu of the gaussian nitrogen density profile
sigma = FWHM/(2*(np.sqrt(2*np.log(2)))) #sigma of the gaussian density profile
r=1

# Gas densities
am_H2 = 100 - am_N2 #Amount of hydrogen
n_tot_gas = 0.5*n_tot*(1/(am_H2/100 + 5*am_N2/100 ))
n_gas_H2 = (am_H2/100)*n_tot_gas #Gas density of H2
n_gas_N2 = (am_N2/100)*n_tot_gas #Gas density of N2


# Initialization of the density arrays
a0_ts = np.zeros_like(ts.t, dtype=object)
gas_dens_H2 = a0_ts.copy()
gas_dens_N2 = a0_ts.copy()
gas_dens_N2_thres = a0_ts.copy()
e_dens_tot = a0_ts.copy()
gas_dens_tot = a0_ts.copy()
gas_dens_sum = a0_ts.copy()

#Drawing the density profiles

for ii, t in enumerate(ts.t):
	if (t*c) < ramp_start:
		e_dens_tot[ii] = 0
		gas_dens_tot[ii] = 0
	elif (t*c) < (ramp_start+ramp_length):
		e_dens_tot[ii] = n_tot*np.sin((0.5*np.pi*(t*c-ramp_start))/(ramp_length))**2
		gas_dens_tot[ii] = n_tot_gas*np.sin((0.5*np.pi*(t*c-ramp_start))/(ramp_length))**2
	elif (t*c) > (ramp_start+ramp_length + plateau):
		e_dens_tot[ii] = n_tot*np.sin((0.5*np.pi*(t*c-ramp_start))/ramp_length)**2
		gas_dens_tot[ii] = n_tot_gas*np.sin((0.5*np.pi*(t*c-ramp_start))/(ramp_length))**2
	elif (t*c) >= (ramp_start + 2*ramp_length + plateau):
		e_dens_tot[ii] = 0
		gas_dens_tot[ii] = 0
	else:
		e_dens_tot[ii] = n_tot
		gas_dens_tot[ii] = n_tot_gas

	gas_dens_N2[ii] = n_gas_N2*np.exp(-((c*t)-mu)**2/(2*sigma**2)) 
	
	gas_dens_N2_thres[ii] = ((inj_thres-1)/7)*n_gas_N2*np.exp(-((c*t)-mu)**2/(2*sigma**2)) 	
	
	gas_dens_H2[ii] = gas_dens_tot[ii]-gas_dens_N2_thres[ii]

	gas_dens_sum[ii] = gas_dens_H2[ii] + gas_dens_N2[ii]	

# Plot the functions

fig, ax = plt.subplots(figsize=(20,10))

plt.xlabel("z [mm]")
plt.ylabel(r"$n_e$")
   
ax.set_xticks(np.arange(0, (2*ramp_start+2*ramp_length+plateau)*1e3, 1) )                                                       
ax.set_xticks(np.arange(0, (2*ramp_start+2*ramp_length+plateau)*1e3, 0.2), minor=True)                                           
ax.set_yticks(np.arange(0, n_tot*10, 1e24))                                                       
ax.set_yticks(np.arange(0, n_tot*10, 0.2e24), minor=True)   

ax.grid(which='major', alpha=0.8, color='grey', linestyle='--', linewidth=0.6, animated='True')
ax.grid(which='minor', alpha=0.6, color='grey', linestyle='--', linewidth=0.5, animated='True')

# Plot the densities
ax.plot(ts.t*c*1e3, e_dens_tot, linestyle='--' ,color='grey', label="$n_e (tot)$")
ax.plot(ts.t*c*1e3, gas_dens_tot, linestyle='--' ,color='red', label="$n_{gas} (tot)$")
ax.plot(ts.t*c*1e3, gas_dens_H2, linestyle='--' ,color='#31B5D6', label="$n_{gas} (H2)$")
ax.plot(ts.t*c*1e3, gas_dens_N2, linestyle='--', color='#7BC618', label="$n_{gas} (N2)$")
ax.plot(ts.t*c*1e3, gas_dens_N2_thres, linestyle='--', color='green', label="$n_{gas,thres} (N2)$")
ax.plot(ts.t*c*1e3, gas_dens_sum, linestyle='--', color='black', label="$n_{gas,real} (tot)$")

legend = ax.legend(loc='upper right', shadow=True)
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
title("Density Profile Analytics\n\n initial gas density: "+str(n_tot)+" [nitrogen: "+str(am_N2)+"%, hydrogen: "+str(100-am_N2)+"%]", bbox={'facecolor': '0.85', 'pad': 10})

show()


