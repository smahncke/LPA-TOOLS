#IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from pylab import title, show
from scipy.constants import e, c, m_e, h, alpha, hbar

#Initial parameters
lambda_0 = 0.8e-6 #laser wavelength in m
initial_energy = 1e3 #initial energy of the ionized electrons (in eV)
E_L = 8.e12 #Field of the laser in V/m
element = 'Ar' # Element that is used fo ionization injection (string, use the shortcut, f.e. 'N' for nitrogen)



#---------------------------------------
# Don't change anything below this line
#---------------------------------------

#Ionization energies (in eV)

if element == 'N':
	U_G = [14.534,29.600,47.449,77.473,97.890,552.070,667.045]
elif element == 'Ar':
	U_G = [15.760,27.629,40.742,59.812,75.017,91.009,124.319,143.462,422.448,478.684,\
			538.963,618.260,686.104,755.742,854.772,918.026,4120.887,4426.228]
elif element == 'Kr':
	U_G = [14.000,24.360,26.949,52.547,64.673,78.458,111.001,125.802,230.854,\
			268.227,307.819,350.312,390.733,446.700,492.303,541.015,591.800,640.512,785.612,833.287,\
			884.072,936.930,998.079,1050.937,1151.471,1205.261,2927,907,3069.897,3227.434,3380.826]
else:
	print("Choose an available element for ionization injection")

#Physical constants
lambda_c = h/(m_e * c) #compton wavelength
U_H = 13.6 #ionization potential of hydrogen in eV

#Other parameters
k_0 = 1/lambda_0 #wave number of the laser
E_k = (1/e)*k_0*(c**2)*m_e #laser energy
omega_0 = c*k_0 #laser frequency
energy = initial_energy*e 

#ionization probability
def ion_prob(U_i,amplitude):
	""" Inputs: 
        - U_i (float): Ionization potential
		- amplitude (1D array): Laser amplitude 
        Output: Probability, that the electron got ionized at rest
				in dependence of the laser amplitude
	"""
	#Keldysh parameter
	gamma_k = (alpha/amplitude)*np.sqrt(U_i/U_H) 
    
    #Argument of the exponential function
	laser_ion = lambda_0/lambda_c * amplitude**3 *gamma_k**3 * (E_k/(E_L))
	tunnel_ion = (gamma_k**3 * energy)/(hbar*omega_0)
    
    #Probability
	prob = np.exp(-2/3 * (laser_ion + tunnel_ion))
    
	return(prob)

def energy_spread(U_i):
	"""	Input: U_i (float): Potential of the ionization level (in eV)
		Output: delta_E (float): Energy spread of the corresponding ionization level
	"""
	#Keldysh parameter	
	gamma_k = (alpha/amplitude)*np.sqrt(U_i/U_H) 

	#Energy spread
	delta_E = 3*hbar*omega_0/(2*e*1e6*gamma_k**3)

	return(delta_E)

#The possible energy range
a0 = np.linspace(0.01,2.5,1000)


### Plot the functions

# Set the size of the plot
fig, ax = plt.subplots(figsize=(20,10))

# Set the labels of the axis
plt.xlabel("Laser amplitude a")
plt.ylabel(r"Probability $w(a)$")

ax.grid(which='major', alpha=0.8, color='grey', linestyle='--', linewidth=0.6, animated='True')
ax.grid(which='minor', alpha=0.6, color='grey', linestyle='--', linewidth=0.5, animated='True')

#Make the plots
for i in range(0,len(U_G)): 
	ax.plot(a0, ion_prob(U_G[i],a0), label=""+str(element)+"$^{"+str(i+1)+"+}$")

# Make the legend 
legend = ax.legend(loc='upper right', shadow=True)
for label in legend.get_texts():
    label.set_fontsize('medium')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
title("Ionization probability [ $\lambda_0 =$ "+str(lambda_0*1.e9)+" nm, $\epsilon =$"+str(energy/e)+" eV , $|E_L| =$"+str(E_L/1.e12)+" TV/m]", bbox={'facecolor': '0.85', 'pad': 10})

# Opens the plot
show()



