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

### Global settings

ts = OpenPMDTimeSeries('./diags/hdf5/') #The simulation files which should be imported

### THE GAS

n_tot = 1e18*1e6  # Total gas density
am_N2 = 10.   #Amount of total N2 molecules (including pre-ionized and non pre-ionized ones) in percent
#pre_ion_level = int(input("Pre-Ion-Level"))
pre_ion_level = 0 #The pre-ionization level of the N2 (0 = no pre-ionization)
max_ion_level = 7 # The highest possible ionization level (7 for nitrogen)
E_ion = [5,10,15,20,25,30,35] # The ionization energy levels of N2
#E_ion = [14.53414, 29.6013, 47.44924, 77.4735, 97.8902, 552.0718, 667.046] #Ionization energies of N2 in eV

# The dope gas density profile (Gauss)

FWHM = 300.e-6 #FWHM of the gauss function in m
ramp_start = 30.e-6 # The starting point of the density profile upramp
ramp_length = 500.e-6 # The length of the ramp
plateau = 4000.e-6 # The length of the plateau



#################################################################################################################
#														#
#				DON'T CHANGE ANYTHING BELOW THIS LINE						#
#														#
#################################################################################################################

print("\n--------------------- Plasma Density Plotter [PDP] ---------------------")
time.sleep(1)
print("\nWelcome!\n")
print("> Current version: 0.1")
print("> visit 'https://github.com/smahncke/pdp' for further information\n------------------------------------------------------------------------\n")
time.sleep(1.5)
print("The current parameters are: \n \n - Total gas density: "+str(n_tot)+"\n - Gas ratio: "+str(am_N2)+"% nitrogen, "+str(100.-am_N2)+"% hydrogen\n - pre-ionization-level: "+str(pre_ion_level)+"\n")


main_menu = input("Press 'ENTER' to start or type 'e' to edit the values\n")

if main_menu == "e":
	print("------------------------------------------------------------------------\n> Parameter menu:\n")
	while 1:
		sub_menu_1 = input("Type in one of the given button to change the corresponding parameter:\n\n n: Total initial gas density\n a: Amount of nitrogen (in percent)\n p: Pre-ion level\n\nType in 'e' to go back ")
		if sub_menu_1 == "n":
			n_tot = input("\nType in a new value for the total initial gas density and press 'ENTER': ")
		elif sub_menu_1 == "a":
			am_N2 = int(input("\nType in a new value for the inital amount of nitrogen gas in percent and press 'ENTER': "))
		elif sub_menu_1 == "p":
			pre_ion_level = int(input("\nType in a new value for the pre-ionization level and press 'ENTER': "))
		elif sub_menu_1 =="e":
			break

	print("------------------------------------------------------------------------\n> Parameters edited!\n\nThe current values are:\n\n- Total initial gas density: "+str(n_tot)+"\n- Initial amount of nitrogen: "+str(am_N2)+"%\n- Pre-ionization level: "+str(pre_ion_level)+"\n")
	not input("Press 'ENTER' to start\n")
print("-----------------------\nOperation started\n-----------------------\n")
	
time.sleep(1.5)
### Calculations for the gas densities profile

# Parameters of the gaussian profile
mu = ramp_start + ramp_length # mu of the gaussian dope gas density profile
sigma = FWHM/(2*(np.sqrt(2*np.log(2)))) # sigma of the gaussian density profile
r=1

# Gas densities
am_H2 = 100. - am_N2 #Amount of (already pre-ionized) H2
n_gas_H2 = (am_H2/100)*n_tot # gas density of H2
n_gas_N2 = (am_N2/100)*n_tot # gas density of N2

# Plasma densities
n_plas_H2 = n_gas_H2*2. #Plasma density of the H2 (without the pre-ionized electrons)
n_plas_N2 = n_gas_N2*2 #Plasma density of the N2 (without pre-ionization)

# Electron densities
n_e_N2_init = max_ion_level*n_plas_N2 #Electron density of the initial N2
n_e_N2_pre_ion = n_plas_N2*pre_ion_level #Electron density of the pre-ionized N2 (goes into the H2 plasma density)
n_e_N2_rest = n_e_N2_init - n_e_N2_pre_ion #Resulting electron density without the already pre-ionized electrons
n_e_N2_rest_ionized = n_e_N2_rest # Initializes the ionized electrons
n_e_tot = n_e_N2_init + n_plas_H2 # Total electron density

#if pre_ion_level > max_ion_level: 
#    raise ValueError("ERROR: pre_ion_level is higher than max_ion_level (Currently at "+str(pre_ion_level)+", the highest allowed value is "+str(max_ion_level)+")"

### DENSITY PROFILES
print("Creating the arrays..")
# Initialize the arrays
a0_ts = np.zeros_like(ts.t , dtype=object) 
delta_n = a0_ts.copy() # The amount of electrons, which goes into the H2 due to laser ionization 
n_e_N2_rest_ionized = a0_ts.copy() 

e_dens_N2_ionized = a0_ts.copy()
e_dens_N2_rest = a0_ts.copy()
e_dens_N2_init = a0_ts.copy()

e_dens_tot = a0_ts.copy()
e_dens_H2 = a0_ts.copy()
e_dens_H2_init = a0_ts.copy()
time.sleep(1.5)
print("Arrays successfully created...\n")
time.sleep(0.5)
print("Creating the plot...")
time.sleep(1)


### IONIZATION & DENSITY ARRAYS

for ii, t in enumerate(ts.t): # A for loop over the whole plasma
    sys.stdout.write("\r" + "Status: " + str(int((ii/(len(ts.t)))*100)+1) + " %")
	
    a0_ts[ii] = ts.get_a0(t=t, pol='x')*10 # Gives the laser amplitude for each time step
    
    # Total density profile
    if (t*c) < ramp_start:
        e_dens_tot[ii] = 0
    elif (t*c) < (ramp_start+ramp_length):
        e_dens_tot[ii] = n_e_tot*np.sin((0.5*np.pi*(t*c-ramp_start))/(ramp_length))**2
    elif (t*c) > (ramp_start+ramp_length + plateau):
        e_dens_tot[ii] = n_e_tot*np.sin((0.5*np.pi*(t*c-ramp_start))/ramp_length)**2
    elif (t*c) >= (ramp_start + 2*ramp_length + plateau):
        e_dens_tot[ii] = 0
    else:
        e_dens_tot[ii] = n_e_tot
    
    for i in range(0,int(max_ion_level)-1): # Compares the laser amplitude to each ionization energy level of N2
        if E_ion[i] <= a0_ts[ii] < E_ion[i+1]:
            delta_n[ii] = n_plas_N2*((i+1)-pre_ion_level)
        elif a0_ts[ii] >= E_ion[6]:
            delta_n[ii] = n_plas_N2*(7-pre_ion_level)

        if delta_n[ii] < 0: # Makes sure that the N2 density doesn't get greater by the laser
            delta_n[ii] = 0
    
        n_e_N2_rest_ionized[ii] = n_e_N2_rest - delta_n[ii] # The resulting N2 electron density AFTER laser ionization

        if n_e_N2_rest_ionized[ii] <= 1e10: # Compensates rounding errors of the density
            n_e_N2_rest_ionized[ii] = 0  
            
        # The initial electron density of N2 (with pre-ionization, but without laser-ionization)
        e_dens_N2_init[ii] = n_e_N2_rest*np.exp(-((c*t)-mu)**2/(2*sigma**2)) 
        
        # The initial electron density of H2 (with pre-ionization, but without laser-ionization)
        e_dens_H2_init[ii] = e_dens_tot[ii] - e_dens_N2_init[ii]
        
        # The electron density, which goes from the N2 to the H2 due to laser ionization:
        e_dens_N2_ionized[ii] = delta_n[ii]*np.exp(-((c*t)-mu)**2/(2*sigma**2)) 
        
        # The rest of the N2 electron density
        e_dens_N2_rest[ii] = n_e_N2_rest_ionized[ii]*np.exp(-((c*t)-mu)**2/(2*sigma**2))
        
        

        
    # H2 density profile
    #e_dens_H2[ii] = e_dens_tot[ii] - e_dens_N2_ionized[ii]
    e_dens_H2[ii] = e_dens_tot[ii] - e_dens_N2_rest[ii]
    sys.stdout.flush()
print()
time.sleep(0.5)

### Plot the functions (plasma densities)

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
ax.plot(ts.t*c*1e3, e_dens_H2_init, linestyle='--' ,color='#31B5D6', label="$n_e (H2\,init)$")
ax.plot(ts.t*c*1e3, e_dens_N2_init, linestyle='--', color='#7BC618', label="$n_e (N2\,init)$")
ax.plot(ts.t*c*1e3, e_dens_H2,color='#005263', label="$n_e (H2)$")
ax.plot(ts.t*c*1e3, e_dens_N2_rest,color='#218429', label="$n_e (N2)$")
ax.plot(ts.t*c*1e3, e_dens_N2_ionized,linestyle='--', color='#FF8429', label="$\Delta n_e$")

#Plot the laser amplitude
ax_las = ax.twinx()
plt.ylabel(r"$a_0$", color="red")
ax_las.plot(ts.t*c*1e3, a0_ts,color='red')

legend = ax.legend(loc='upper right', shadow=True)
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
title("Density Profile Analytics\n\n initial gas density: "+str(n_tot)+" [nitrogen: "+str(am_N2)+"%, hydrogen: "+str(am_H2)+"%, pre-ionization level: "+str(pre_ion_level)+"]", bbox={'facecolor': '0.85', 'pad': 10})

print("Plot successfully created...\n")
time.sleep(1)
save_plot = input("Save the plot as a .png file? y/[n]")
if save_plot == "y":
	if not os.path.exists(os.path.dirname("plots/")):
		print("Folder 'plots/' not found, so it will be created!\nThe plot will be saved in this folder..")
		os.makedirs(os.path.dirname("plots/"))
		time.sleep(1)
		print("Done!\n")
	else:
		print("The folder 'plots/' does already exist!\nThe plot will be saved in this folder..\n")
	time.sleep(0.5)
	print("Saving as .png file...")
	fig.savefig("plots/dpa_plot_N2_"+str(int(am_N2))+"_preionlvl_"+str(pre_ion_level)+".png")
	time.sleep(1)
	print("Done!\n")

	time.sleep(0.5)

show_plot = input("Open the plot now? y/[n]")
if show_plot == "y":
	print("Opening the plot..")
	time.sleep(1)
	show()

print("-----------------------\nFinished\n-----------------------")

