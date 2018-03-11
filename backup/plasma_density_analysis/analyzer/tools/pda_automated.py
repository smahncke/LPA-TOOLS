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

### IMPORT THE DIAGNOSTICS

ts = OpenPMDTimeSeries('./diags/hdf5/') #The simulation files which should be imported

### INITIAL GLOBAL PARAMETERS

n_tot = 1e18*1e6  # Total gas density

am_N2 = 1   #Amount of total initial N2 molecules (including pre-ionized and non pre-ionized ones) in percent
max_am_N2 = 2 #The amount of N2 at which the automated script should stop (in percent)
am_N2_step = 1 #The loop step of the simulation (in percent)

pre_ion_level = 0 #The pre-ionization level of the N2 (0 = no pre-ionization)
max_pre_ion_level = 2 #The highest value for pre-ionization the plots will be created for
max_ion_level = 7 # The highest possible ionization level (7 for nitrogen)
E_ion = [5,10,15,20,25,30,35] # The ionization energy levels of N2
#E_ion = [14.53414, 29.6013, 47.44924, 77.4735, 97.8902, 552.0718, 667.046] #Ionization energies of N2 in eV

### THE DENSITY PROFILES

FWHM = 300.e-6 #FWHM of the gauss function in m
ramp_start = 30.e-6 # The starting point of the density profile upramp
ramp_length = 500.e-6 # The length of the ramp
plateau = 4000.e-6 # The length of the plateau



#####################################################################
#																	#
#				DON'T CHANGE ANYTHING BELOW THIS LINE				#
#																	#
#####################################################################

### START SCREEN


while 1:
    screen_total_gas_dens = " - Total gas density: "+str(n_tot)+""

    if am_N2 == max_am_N2:
        screen_gas_ratio = " - Constant gas ratio: " +str(am_N2)+"%"
    else:
        screen_gas_ratio = " - Varying gas ratio in range of " +str(am_N2)+"% to "+str(max_am_N2)+"% nitrogen in "+str(am_N2_step)+"% steps"

    if pre_ion_level == max_pre_ion_level:
        if pre_ion_level == 0:
            screen_pre_ion_level = " - Constant pre-ionization level: 0 (no pre-ionization)"
        else:
            screen_pre_ion_level = " - Constant pre-ionization level:" +str(pre_ion_level)+""
    else:
        screen_pre_ion_level = " - Varying pre-ionization level in range of " +str(pre_ion_level)+ " to " + str(max_pre_ion_level)+""
    
    print("\n---------------- Plasma Density Plotter ----------------\n")
    print("Current simulaton settings:")
    print(screen_total_gas_dens)
    print(screen_gas_ratio)
    print(screen_pre_ion_level)
    print("----------------------------------------------------------\n")

    main_menu = input("> Press 'ENTER' to start or type 'e' to edit the values\n")

    ### EDIT MENU

    if main_menu == "e":
        print("----------------------------------------\nParameter menu:\n----------------------------------------")
        while 1:
            sub_menu_1 = input("> Type in one of the given button to change the corresponding parameter:\n\nn: Total initial gas density\na: Amount of nitrogen (in percent)\np: Pre-ion level\n\nType in 'e' to go back ")
            if sub_menu_1 == "n":
                n_tot = input("\n> Type in a new value for the total initial gas density and press 'ENTER': ")
            elif sub_menu_1 == "a":
                am_N2 = int(input("\n> Initial amount of nitrogen (in %): "))
                max_am_N2 = int(input("> Final amount of nitrogen (in %): "))
                am_N2_step = int(input("> Stepwidth (in %): "))
                if max_pre_ion_level < pre_ion_level:
                    max_pre_ion_level = pre_ion_level
                    print("\nThe highest pre-ionization level must be at least as high as the initial pre-ionization level")
            elif sub_menu_1 == "p":
                pre_ion_level = int(input("\n> Initial pre-ionization level (0 = no pre-ionization): "))
                max_pre_ion_level = int(input("> Final ionization level: "))
            elif sub_menu_1 =="e":
                print("Saved changes!")
                break
    else:
        break

init_pre_ion_level = pre_ion_level
print("----------------------------------------\nOperation started\n----------------------------------------\n")
	
time.sleep(0.5)


###### MAIN SCRIPT ######

mu = ramp_start + ramp_length # mu of the gaussian dope gas density profile
sigma = FWHM/(2*(np.sqrt(2*np.log(2)))) # sigma of the gaussian density profile
r=1

## Initialize the loop parameters
number_of_plots = int((max_pre_ion_level+1)-pre_ion_level)*(((max_am_N2+1)-am_N2)/am_N2_step) # The total number of plots
current_plot = 0 # Gives the number of the current plot

print("Total number of plots: "+str(int(number_of_plots))+"...\n") #Prints out the total number of plots

time.sleep(1)

## LOOP 1: Varies the amount of nitrogen 
for y in range(int(am_N2),max_am_N2+1):
    
    #print("\nCurrent amount of nitrogen: "+str(am_N2)+"%")
    if not os.path.exists(os.path.dirname("plots/N2_"+str(am_N2)+"/")): #Creates the directory in which the plots will be saved (if it doesn't already exist)
        os.makedirs(os.path.dirname("plots/N2_"+str(am_N2)+"/"))

	# Gas densities
    am_H2 = 100. - am_N2 #Amount of (already pre-ionized) H2
    n_gas_H2 = (am_H2/100)*n_tot # gas density of H2 
    n_gas_N2 = (am_N2/100)*n_tot # gas density of N2

	# Plasma densities
    n_plas_H2 = n_gas_H2*2. #Plasma density of the H2 (without the pre-ionized electrons)
    n_plas_N2 = n_gas_N2*2 #Plasma density of the N2 (without pre-ionization)
	
	## LOOP 2: Varies the pre-ionization level
    for x in range(pre_ion_level,max_pre_ion_level+1):

        n_e_N2_init = max_ion_level*n_plas_N2 #Electron density of the initial N2
        n_e_N2_pre_ion = n_plas_N2*pre_ion_level #Electron density of the pre-ionized N2 (goes into the H2 plasma density)
        n_e_N2_rest = n_e_N2_init - n_e_N2_pre_ion #Resulting electron density without the already pre-ionized electrons
        n_e_N2_rest_ionized = n_e_N2_rest # Initializes the ionized electrons
        n_e_tot = n_e_N2_init + n_plas_H2 # Total electron density

	# Initialize the arrays
        a0_ts = np.zeros_like(ts.t , dtype=object) 
        delta_n = a0_ts.copy() # The amount of electrons, which goes into the H2 due to laser ionization 
        n_e_N2_rest_ionized = a0_ts.copy() 

        e_dens_N2_ionized = a0_ts.copy() # electron density array of the electrons, which go from the nitrogen to the hydrogen due to laser ionization
        e_dens_N2_rest = a0_ts.copy() # electron density array of resulting nitrogen gas after laser ionization
        e_dens_N2_init = a0_ts.copy() # electron density array of the initial nitrogen gas

        e_dens_tot = a0_ts.copy() # Total electron density
        e_dens_H2 = a0_ts.copy() # electron density array of the resulting hydrogen gas after laser ionization
        e_dens_H2_init = a0_ts.copy() # electron density array of the initial hydrogen gas
		
		### LOOP 3: Creates the plot over the whole plasma length 
        for ii, t in enumerate(ts.t): 
            sys.stdout.write("\r" + "Total Status: " + str(int((current_plot/number_of_plots)*100)) + "%        [ Current Plot: "+str(int(current_plot+1))+" of "+str(int(number_of_plots))+" ("+str(int((ii/(len(ts.t)))*100)+1) + " %" + ") ]")
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
			
			### LOOP 4: Checks for every timestep if the laser energy is higher than the ionization levels calculates the resulting density profiles
            for i in range(0,int(max_ion_level)-1): 
                if E_ion[i] <= a0_ts[ii] < E_ion[i+1]:
                    delta_n[ii] = n_plas_N2*((i+1)-pre_ion_level)
                elif a0_ts[ii] >= E_ion[6]:
                    delta_n[ii] = n_plas_N2*(7-pre_ion_level)

				# Makes sure that the N2 density doesn't get greater by the laser
                if delta_n[ii] < 0: 
                    delta_n[ii] = 0
				
				# The resulting N2 electron density AFTER laser ionization
                n_e_N2_rest_ionized[ii] = n_e_N2_rest - delta_n[ii] 

				# Compensates rounding errors of the density
                if n_e_N2_rest_ionized[ii] <= 1e10: 
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
            e_dens_H2[ii] = e_dens_tot[ii] - e_dens_N2_rest[ii]
            
        #sys.stdout.flush()
		### MAKING THE PLOTS
        fig, ax = plt.subplots(figsize=(20,10))
        plt.xlabel("z [mm]")
        plt.ylabel(r"$n_e$")

		#ax.plot(ts.t*c*1e3, a0_ts,color='red')

		#ax.twinx()

        ax.plot(ts.t*c*1e3, e_dens_tot, linestyle='--' ,color='grey', label="$n_e (tot)$")
        ax.plot(ts.t*c*1e3, e_dens_H2_init, linestyle='--' ,color='#31B5D6', label="$n_e (H2\,init)$")
        ax.plot(ts.t*c*1e3, e_dens_N2_init, linestyle='--', color='#7BC618', label="$n_e (N2\,init)$")
        ax.plot(ts.t*c*1e3, e_dens_H2,color='#005263', label="$n_e (H2)$")
        ax.plot(ts.t*c*1e3, e_dens_N2_rest,color='#218429', label="$n_e (N2)$")
        ax.plot(ts.t*c*1e3, e_dens_N2_ionized,linestyle='--', color='#FF8429', label="$\Delta n_e$")

        legend = ax.legend(loc='upper right', shadow=True) 
        for label in legend.get_texts():
            label.set_fontsize('large')

        for label in legend.get_lines():
            label.set_linewidth(1.5)  # the legend line width
        title("Density Profile Analytics\n\n initial gas density: "+str(n_tot)+" [nitrogen: "+str(am_N2)+"%, hydrogen: "+str(am_H2)+"%, pre-ionization level: "+str(pre_ion_level)+"]", bbox={'facecolor': '0.85', 'pad': 10})

		# Saves the plots as a png file
        fig.savefig("plots/N2_"+str(am_N2)+"/dpa_plot_N2_"+str(int(am_N2))+"_preionlvl_"+str(pre_ion_level)+".png")
		
		# Makes sure, that the CPU doesn't get out of memory
        plt.close()
		
        pre_ion_level += 1
        current_plot += 1
    pre_ion_level = init_pre_ion_level
    am_N2 += am_N2_step
    #print()
    sys.stdout.flush()
print()
print("\n----------------------------------------\nFinished\n----------------------------------------")
