#-------
#IMPORTS
#-------

import numpy as np
import time
from scipy.integrate import odeint
from pylab import *
import scipy.constants as const
from numba import double
from numba.decorators import jit, autojit
from scipy.constants import m_e,c,e, epsilon_0, h, alpha, hbar
import matplotlib.pyplot as plt
from pylab import title, show

#-----
#TIMER
#-----

start_timer = time.time()

#------------------
#INITIAL PARAMETERS
#------------------

n0 = 1.e24 #Plasma density
lambda_0 = 800.e-9 #laser wavelength
RMS = 1 #RMS of the gaussian laser profile
E_max = 2 

a0 = 1.2 #Max laser amplitude

w0 = 35.e-6 #Laser waist
ctau = 14.e-6 #Laser duration
z0 = 0.e-6 #Focus position
zf = 500.e-6 #Laser focus

#GAS
energy = 0e3 #initial energy of the ionized electrons (in eV)
element = 'Ar' # Element that is used fo ionization injection (string, use the shortcut, f.e. 'N' for nitrogen)

#PLOT SETTINGS
ion_niveau = 8 #The niveau that should get ionized 	
iterations = 1.2e7 #Resolution of the plots (the higher the better)

#------------------------------
#OTHER PARAMETERS AND CONSTANTS
#------------------------------

#Physical constants
lambda_c = h/(m_e * c) #compton wavelength
U_H = 13.6 #ionization potential of hydrogen in eV
epsilon_0 = 8.854e-12 #Vacuum permittivity

#Other parameters
omega_0 = 2*np.pi*c/lambda_0 # Omega of the laser
k_0 = 2*np.pi/lambda_0 #wave number of the laser
E_k = (1/e)*k_0*(c**2)*m_e #laser energy

# Density distribution
def n(zz):  
    return n0 * ones(zz.size)

# Plasma frequency
def omega_p(zz): 
    return sqrt(n(zz)*e**2/(epsilon_0*m_e))

# Initial electric field of the plasma in V/m
def E_0(zz):
    return c*m_e*omega_p(zz)/e 

# Normalized plasma phase velocity
def beta_p(zz):
    return 1-omega_p(zz)**2/(2*omega_0**2)
    
# Plasma wave phase velocity
def gamma_p(zz):
    return np.sqrt(1/(1-beta_p(zz)**2))

# Plasma wavenumber
def k_p(zz): 
    return omega_p(zz)/c

#----------------------
#INITIALIZE OUTPUT LIST
#----------------------

output_list = ['Ionization induced trapping - outputs\n-----------------------------\n']
output_list.append('Element: '+str(element)+'')
output_list.append('Ionization level: '+str(element)+str(ion_niveau)+'+ -> '+str(element)+str(ion_niveau+1)+'+')
output_list.append('Max laser amplitude a0: '+str(a0)+'\n')

#----------
#INFO TEXTS
#----------

print('\n------------------------------\nCALCULATIONS STARTED\n')

#--------------
#LASER PROFILES
#--------------

#Envelope
def a(psi,a0,RMS):
    return a0*np.exp(-psi**2/(4*RMS**2))

#Oscillating field
def gaussian_profile(max_a, z, r, t, w0, ctau, z0, zf, k_0):

    # Calculate the Rayleigh length

    zr = 0.5*k_0*w0**2
    inv_zr = 1./zr
    inv_ctau2 = 1./ctau**2
    
    
    # Diffraction and stretch_factor
    diffract_factor = 1. - 1j*(z-zf)*inv_zr
    stretch_factor = 1.
    # Calculate the argument of the complex exponential
    exp_argument = 1j*k_0*( c*t + z0 - z )         - r**2 / (w0**2 * diffract_factor)         - 1./stretch_factor * inv_ctau2 * ( c*t  + z0 - z )**2
    # Get the transverse profile
    profile_Eperp = np.exp(exp_argument)         / ( diffract_factor * stretch_factor**0.5 )
    
    return(max_a*profile_Eperp.real )


#------------------
#IONIZATION PROCESS
#------------------

#IONIZATION PROBABILITY
def ion_prob(psi):
    
    E_gauss = gaussian_profile(a0, psi,0,0,w0,ctau,z0,zf,k_0) 
    
    amplitude = a(psi,a0, RMS)

	#Keldysh parameter
    gamma_k = (alpha/amplitude)*np.sqrt(U_i/U_H) 
    
	#Calculate the electrical field of the laser (in V/m)
    E_L = 1e2*np.sqrt(2.8e18/(lambda_0*1e6*c*epsilon_0))*E_gauss #Electrical field of the laser in V/m
	
    #Argument of the exponential function
    laser_ion = lambda_0/lambda_c * amplitude**3 *gamma_k**3 * (E_k/(np.sqrt(E_L**2)))
    tunnel_ion = (gamma_k**3 * energy)/(hbar*omega_0)
    
    #Probability
    prob = np.exp(-2/3 * (laser_ion + tunnel_ion))
    
    return(prob)

#IONIZATION DEGREE
def degree(psi,info):
	n = np.zeros_like(psi+1)
	print(info)
	#start_timer_2 = time.time()
	for i,val in enumerate(psi):
		sys.stdout.write("\r" + " > Current status:" + str("% .1f" % (i/len(psi)*100)) +"%")
		prob = ion_prob(val)
		if i+1 < len(psi):
			n[i+1]=  n[i]+prob*(n0-n[i])
		sys.stdout.flush()
	print()
	#end_timer_2 = time.time()
	#print("This process took "+str(int((end_timer_2-start_timer_2)/60))+" min, "+str(int(((end_timer_2-start_timer_2)/60)-int((end_timer_2-start_timer_2)/60))*60)+" sec!\n")
	if info == "Calculating ionization degree...":
		output_list.append('- Max. ionization degree: '+ str(n[len(psi)-1]/n0 *100) +'%')
		output_list.append('- Max. ionization probability: '+ str(ion_prob(psi[int(len(psi)/2)])*100)+'%')
		print('- Max. ionization degree: '+ str(n[len(psi)-1]/n0 *100) +'%')
		print('- Max. ionization probability: '+ str(ion_prob(psi[int(len(psi)/2)])*100)+'%\n')
	return n/n0


#IONIZATION ENERGIES
if element == 'N':
	U_G = [14.534,29.600,47.449,77.473,97.890,552.070,667.045]
elif element == 'Ar':
	U_G = [15.760,27.629,40.742,59.812,75.017,91.009,124.319,143.462,422.448,478.684,	538.963,618.260,686.104,755.742,854.772,918.026,4120.887,4426.228]
elif element == 'Kr':
	U_G = [14.000,24.360,26.949,52.547,64.673,78.458,111.001,125.802,230.854,268.227,307.819,350.312,390.733,446.700,492.303,541.015,591.800,640.512,785.612,833.287,884.072,936.930,998.079,1050.937,1151.471,1205.261,2927,907,3069.897,3227.434,3380.826]
else:
	print("Choose an available element for ionization injection")
    
U_i = U_G[ion_niveau]

          
#----------------------------------
#QUASI STATIC PLASMA WAVE EEQUATION          
#----------------------------------
          
# Potential (PSI DEPENDENCE)
def potentialeq(phi,psi,z,zi):
    return array([phi[1],gamma_p(z)[zi]**2 * (beta_p(z)[zi]*(np.sqrt(1/(1 - (1+ a(psi,a0,RMS)**2)/(gamma_p(z)[zi]**2 * (1+phi[0])**2))))-1)])

# Phi min
def get_phi_min(E_max,z):
    return (-1+1/(2*gamma_p(z)**2))*(1-1/E_max**2) + E_max**2/(4*gamma_p(z)**2)

#-----
#PLOTS
#-----

fig, ax = plt.subplots(figsize=(20,10))


psi = np.linspace(-2*np.pi,2*np.pi,iterations)
z = np.linspace(0, 1000, iterations)
    
# Initial conditions
phiinit = array([0,0]) 

start_timer_1 = time.time()

print("Calculating the wakefield potential, this may take a moment...")    
# Solve the quasi static plasma wave equation
phi = list(map(lambda x: odeint(potentialeq, phiinit, psi, args = (z,x,)),range(10)))
end_timer_1 = time.time()
print("Finished - This process took "+str(int((end_timer_1-start_timer_1)/60))+" min, "+str(int((((end_timer_1-start_timer_1)/60)-int((end_timer_1-start_timer_1)/60))*60))+" sec!\n")

# Plot of the solution
#ax.plot(psi/np.pi,phi[3][:,0], linestyle='--', label="$a_0 =$ "+str("%.1f" % a0)+"")

# Plot of the trapping condition
    
def left_side(zz): # Left side of the trapping equation
    return 1 + get_phi_min(E_max,zz) - phi[3][:,0]

def right_side(zz): # Right side of the trapping equation
    return np.sqrt(1+a(psi,a0,RMS)**2)/gamma_p(zz)

def condition_fullfilled(zz): # Checks if the trapping condition is fullfilled (right_side() >= left_side())
    diff = np.zeros_like(psi, dtype=object)
    diff = np.where(right_side(zz)-left_side(zz) < 0,0,right_side(zz)-left_side(zz))
    return diff # Just =/= 0 if condition is fullfilled

def trapping(zz,info):
	condition = condition_fullfilled(zz)
	grad_deg = np.gradient(degree(zz/1.4e5,info))
	trap = np.zeros_like(psi+1)
	for ii,val in enumerate(psi):
		if ii+1 < len(psi):
			if condition[ii] != 0:
				trap[ii+1]=  trap[ii]+grad_deg[ii]
			else:
				trap[ii+1] = trap[ii]
	output_list.append('- Max. trapping: '+ str(trap[len(psi)-1] *100) +'%')
	return trap
        
    
# Plot the condition 
ax.plot(psi/np.pi,condition_fullfilled(z),label="Trapping condition")
    
# Plot the left and the right side of the trapping condition equation
#ax.plot(psi/np.pi,right_side(z)-left_side(z),label="$a_0 =$ "+str("%.1f" % a0)+"")
#ax.plot(psi/np.pi,a(psi,a0,RMS), linestyle='--')

# Set the labels
ax.set_xlabel(r'$\psi\,[\pi]$',fontdict={'fontsize':20})
ax.set_ylabel('$H_i - H_S$',fontdict={'fontsize':20})

# Plot an example of the laser profile
ax.plot(psi/np.pi,a(psi,a0,RMS), linestyle='--', color="grey", label="Envelope ($a_0 =$ "+str("%.1f" % a0)+")")
ax.plot(psi/(np.pi), gaussian_profile(a0, psi/1.4e5, 0, 0, w0, ctau, z0, zf, k_0),linestyle='--', label="Laser")
ax.plot(psi/np.pi, degree(psi/1.4e5, "Calculating ionization degree..."), label="Inoization degree")
ax.plot(psi/np.pi, ion_prob(psi/1.4e5), label="Ionization probability")
ax.plot(psi/np.pi, trapping(psi, "Calculating trapping condition..."), label="Trapping")

# Make the legend
legend = ax.legend(loc='upper left', shadow=True)
for label in legend.get_texts():
    label.set_fontsize('medium')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
title("Trapping degree [ $a_0 =$"+str(a0)+", level "+str(element)+"$^{"+str(ion_niveau)+"+}\mapsto$ "+str(element)+"$^{"+str(ion_niveau+1)+"+}$]", bbox={'facecolor': '0.85', 'pad': 10})

          
fig.savefig("plt_"+str(element)+str(ion_niveau)+".png")

print("> Saved plot as 'plt_"+str(element)+str(ion_niveau)+".png'!")

print("\n> Exports results as .dat file\n")


with open('output.dat', 'w') as output_dat:
	for output_list_entries in output_list:
		output_dat.write(str(output_list_entries)+ "\n")

end_timer = time.time()

print("Finished - The calculations took "+str(int((end_timer-start_timer)/60))+" min, "+str(int((((end_timer-start_timer)/60)-int((end_timer-start_timer)/60))*60))+" sec!")

plt.show()
