
### IMPORTS CLASSES
import numpy as np
import csv
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
am_N2 = 30 #Amount of nitrogen (in percent)

# Density profiles
FWHM = 300e-6 #Full-width at half-max of the gaussian density profile (in m)
ramp_start = 30e-6 #Start point of the resulting plasma upramp (in m)
ramp_length = 500e-6 #The length of the upramp and the downramp (in m)
plateau = 4000e-6 #The length of the plateau (in m)

# Other settings
inj_thres = 6 #The ionization level of N2, at which the electron injection starts

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
a0_ts = np.zeros_like(ts.t, dtype=object)
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


### Calculating the density profiles
with open('output/N2_'+str(am_N2)+'/H2_density.csv', 'w', newline='\n') as H2_array, open('output/N2_'+str(am_N2)+'/N2_density.csv', 'w', newline='\n') as N2_array, open('output/N2_'+str(am_N2)+'/total_gas_density.csv', 'w', newline='\n') as total_gas_array:
    
    # Initialization of the CSV tables
    H2_array_writer = csv.writer(H2_array, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    N2_array_writer = csv.writer(N2_array, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    total_gas_writer = csv.writer(total_gas_array, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    # Write the header lines into the CSV files
    H2_array_writer.writerow(['H2 density'])
    N2_array_writer.writerow(['N2 density'])
    total_gas_writer.writerow(['Total gas density'])
    
    
    for ii, t in enumerate(ts.t):
        
        # Density of the main plasma
        if (t*c) < ramp_start:
            e_dens_tot[ii] = 0
        elif (t*c) < (ramp_start+ramp_length):
            e_dens_tot[ii] = n_tot*np.sin((0.5*np.pi*(t*c-ramp_start))/(ramp_length))**2
        elif (t*c) > (ramp_start+ramp_length + plateau):
            e_dens_tot[ii] = n_tot*np.sin((0.5*np.pi*(t*c-ramp_start))/ramp_length)**2
        elif (t*c) >= (ramp_start + 2*ramp_length + plateau):
            e_dens_tot[ii] = 0
        else:
            e_dens_tot[ii] = n_tot
		
        #???
        unknown_dens[ii] = 0.5*e_dens_tot[ii]*(1/(1-(am_N2/100 -(inj_thres-1)*(am_N2/100))*np.exp(-((c*t)-mu)**2/(2*sigma**2))))
        
        # Gas density of the nitrogen
        gas_dens_N2[ii] = (unknown_dens[ii]-0.5*e_dens_tot[ii])*(1/(2-inj_thres))*np.exp(-((c*t)-mu)**2/(2*sigma**2)) 
        
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
        
        # Write the data into the CSV files
        H2_array_writer.writerow([str(gas_dens_H2[ii])])
        N2_array_writer.writerow([str(gas_dens_N2[ii])])
        total_gas_writer.writerow([str(gas_dens_tot[ii])])


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
ax.plot(ts.t*c*1e3, e_dens_tot, linestyle='--' ,color='grey', label="$n_e (tot)$")
ax.plot(ts.t*c*1e3, gas_dens_tot, linestyle='--' ,color='red', label="$n_{gas} (tot)$")
ax.plot(ts.t*c*1e3, gas_dens_H2 ,color='#31B5D6', label="$n_{gas} (H2)$")
ax.plot(ts.t*c*1e3, gas_dens_N2, linestyle='--', color='#7BC618', label="$n_{gas} (N2)$")

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


