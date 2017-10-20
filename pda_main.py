
### IMPORTS CLASSES
import numpy as np
import math
import os
import matplotlib.pyplot as plt
from pylab import title, show
from scipy.constants import c

# Gas settings
n_tot = 1e18*1e6 #Total Plasma Density
am_N2 = 100 #Amount of nitrogen (in percent)

# Density profiles
FWHM = 300e-6 #Full-width at half-max of the gaussian density profile (in m)
ramp_start = 30e-6 #Start point of the resulting plasma upramp (in m)
ramp_length = 500e-6 #The length of the upramp and the downramp (in m)
plateau = 4000e-6 #The length of the plateau (in m)

# Other settings
inj_thres = 6 #The ionization level of N2, at which the electron injection starts

iterations = 10000 #Number of iterations (higher values result in a smoother plot, minimum should be 1000)

####################################################################################
# - Dont change anything below this line


#OS
if not os.path.exists(os.path.dirname('output/N2_'+str(am_N2)+'/')): #Creates the output folder
    os.makedirs(os.path.dirname('output/N2_'+str(am_N2)+'/'))

# Gaussian profile of the nitrogen
mu = ramp_start + ramp_length #mu of the gaussian nitrogen density profile
sigma = FWHM/(2*(np.sqrt(2*np.log(2)))) #sigma of the gaussian density profile
r=1

am_H2 = 100 - am_N2 #Amount of hydrogen

# Initialization of the density arrays
z_array = np.linspace(0,100,num=iterations*1.01)
a0_ts = np.zeros_like(z_array, dtype=object)
gas_dens_H2 = a0_ts.copy()
gas_dens_N2 = a0_ts.copy()
unknown_dens = a0_ts.copy()
e_dens_N2 = a0_ts.copy()
e_dens_H2 = a0_ts.copy()
e_dens_tot = a0_ts.copy()
gas_dens_tot = a0_ts.copy()

# Control arrays
e_dens_tot_test = a0_ts.copy()
gas_dens_H2_test = a0_ts.copy()

dividor = 2e4*(iterations/100)
### Calculating the density profiles

for ii,t in enumerate(z_array):
        
    # Density of the main plasma
    if (ii/dividor) < ramp_start:
        e_dens_tot[ii] = 0
    elif (ii/dividor) < (ramp_start+ramp_length):
        e_dens_tot[ii] = n_tot*np.sin((0.5*np.pi*((ii/dividor)-ramp_start))/(ramp_length))**2
    elif (ramp_start + 2*ramp_length + plateau)>(ii/dividor)>(ramp_start+ramp_length + plateau):
        e_dens_tot[ii] = n_tot*np.sin((0.5*np.pi*((ii/dividor)-ramp_start))/ramp_length)**2
    elif (ii/dividor) >= (ramp_start + 2*ramp_length + plateau):
        e_dens_tot[ii] = 0
    else:
        e_dens_tot[ii] = n_tot
		
    #???
    unknown_dens[ii] = 0.5*e_dens_tot[ii]*(1/(1-(am_N2/100 -(inj_thres-1)*(am_N2/100))*np.exp(-((ii/dividor)-mu)**2/(2*sigma**2))))
        
    # Gas density of the nitrogen
    gas_dens_N2[ii] = (unknown_dens[ii]-0.5*e_dens_tot[ii])*(1/(2-inj_thres))*np.exp(-((ii/dividor)-mu)**2/(2*sigma**2)) 
        
    # Electron density of the nitrogen
    e_dens_N2[ii] = 2*gas_dens_N2[ii]*(inj_thres-1)
        
    # Electron density of the hydrogen
    e_dens_H2[ii] = e_dens_tot[ii] - e_dens_N2[ii]
        
    # Gas density of the hydrogen
    gas_dens_H2[ii] = e_dens_H2[ii]/2
        
    # The total gas density
    gas_dens_tot[ii] = gas_dens_H2[ii]+gas_dens_N2[ii]
        
    # Control functions
    e_dens_tot_test[ii] = 2*gas_dens_H2[ii] + 2*(inj_thres-1)*gas_dens_N2[ii]
    gas_dens_H2_test[ii] = (e_dens_tot[ii] - 2*(inj_thres-1)*gas_dens_N2[ii])/2 
    
        

# Save the arrays as a .dat file
np.savetxt('output/N2_'+str(am_N2)+'/N2_density.dat', gas_dens_N2, delimiter=',')
np.savetxt('output/N2_'+str(am_N2)+'/H2_density.dat', gas_dens_H2, delimiter=',')
np.savetxt('output/N2_'+str(am_N2)+'/tot_gas_density.dat', gas_dens_tot, delimiter=',')

### Plot the functions

# Set the size of the plot
fig, ax = plt.subplots(figsize=(20,10))

# Set the labels of the axis
plt.xlabel("z [mm]")
plt.ylabel(r"$n_e$")

# Draw the grid
ax.set_xticks(np.arange(0, (2*ramp_start+2*ramp_length+plateau)*1e3, 1) )                                                       
ax.set_xticks(np.arange(0, (2*ramp_start+2*ramp_length+plateau)*1e3, 0.2), minor=True)                                           
ax.set_yticks(np.arange(0, n_tot*10, 1e24))                                                       
ax.set_yticks(np.arange(0, n_tot*10, 0.2e24), minor=True)   

ax.grid(which='major', alpha=0.8, color='grey', linestyle='--', linewidth=0.6, animated='True')
ax.grid(which='minor', alpha=0.6, color='grey', linestyle='--', linewidth=0.5, animated='True')

# Make the plots
ax.plot(z_array/20, e_dens_tot, linestyle='--' ,color='grey', label="$n_e (tot)$")
ax.plot(z_array/20, gas_dens_tot, linestyle='--' ,color='red', label="$n_{gas} (tot)$")
ax.plot(z_array/20, gas_dens_H2 ,color='#31B5D6', label="$n_{gas} (H2)$")
ax.plot(z_array/20, gas_dens_N2, linestyle='--', color='#7BC618', label="$n_{gas} (N2)$")

# Test plots
#ax.plot(ts.t*c*1e3, e_dens_tot_test, color='red', label="Test")
#ax.plot(ts.t*c*1e3, gas_dens_H2_test,linestyle='--', color='black', label="Test")

# Make the legend 
legend = ax.legend(loc='upper right', shadow=True)
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
title("Plasma Density Analyzer\n\n Plasma density: "+str(n_tot)+" [nitrogen: "+str(am_N2)+"%, hydrogen: "+str(100-am_N2)+"%]", bbox={'facecolor': '0.85', 'pad': 10})

# Saves the plot as a png file
fig.savefig("output/N2_"+str(am_N2)+"/dpa_plot_N2_"+str(am_N2)+".png")

# Opens the plot
show()


