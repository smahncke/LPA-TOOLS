"""
- RESPOSITORY:	LPA-II TOOLS
- CATEGORY: 	TUNNELING IONIZATION PROBABILITY
- DATE:			10/30/2017

- INFO:	Plots a 2D plot that shows the probability of tunneling ionization for an electron that
		will be born at rest in dependence of the laser amplitude.		

"""

#IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
from pylab import title, show
from scipy.constants import e, c, m_e, h, alpha, hbar

#-------------------
# Initial parameters
#-------------------

#LASER
lambda_0 = 800e-9 #laser wavelength in m
FWHM = 300.e-6 #Full width at  half max of the gaussian laser
a0 = 1.2 #The peak amplitude of the laser

#GAS
initial_energy = 0e3 #initial energy of the ionized electrons (in eV)
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
epsilon_0 = 8.854e-12 #Vacuum permittivity

#Other parameters
k_0 = 1/lambda_0 #wave number of the laser
E_k = (1/e)*k_0*(c**2)*m_e #laser energy
omega_0 = c*k_0 #laser frequency
energy = initial_energy*e 

def gauss_laser(max_amplitude,FWHM,t):
    sigma = FWHM/(np.sqrt(2*np.log(2)))
    laser_a = max_amplitude*np.exp(-(c*t)**2/(2*sigma**2))
    return(laser_a)


#ionization probability
def ion_prob(U_i,max_amplitude,FWHM,t):
    
    amplitude = gauss_laser(max_amplitude, FWHM, t) 
        
	#Keldysh parameter
    gamma_k = (alpha/amplitude)*np.sqrt(U_i/U_H) 
    
	#Calculate the electrical field of the laser (in V/m)
    E_L = 1e2*np.sqrt(2.8e18/(lambda_0*1e6*c*epsilon_0))*amplitude #Electrical field of the laser in V/m
	
    #Argument of the exponential function
    laser_ion = lambda_0/lambda_c * amplitude**3 *gamma_k**3 * (E_k/(E_L))
    tunnel_ion = (gamma_k**3 * energy)/(hbar*omega_0)
    
    #Probability
    prob = np.exp(-2/3 * (laser_ion + tunnel_ion))
    
    return(prob)

t = np.linspace(-10.e-12,10.e-12,1000)


### Plot the functions

# Set the size of the plot
fig, ax = plt.subplots(figsize=(20,10))

# Set the labels of the axis
plt.xlabel("z [mm]")
plt.ylabel(r"Probability $w(a)$")

ax.grid(which='major', alpha=0.8, color='grey', linestyle='--', linewidth=0.6, animated='True')
ax.grid(which='minor', alpha=0.6, color='grey', linestyle='--', linewidth=0.5, animated='True')

#Make the plots
#for i in range(0,len(U_G)): 
    #ax.plot(a0, ion_prob(U_G[i],a0), label=""+str(element)+"$^{"+str(i+1)+"+}$")
ax.plot(c*t, gauss_laser(a0,FWHM,t), label="Laser")
ax.plot(c*t, ion_prob(U_G[5],a0,FWHM,t), label="Probability")
# Make the legend 
legend = ax.legend(loc='upper right', shadow=True)
for label in legend.get_texts():
    label.set_fontsize('medium')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
title("Ionization probability [ $\lambda_0 =$ "+str(lambda_0*1.e9)+" nm, $a_0 =$"+str(a0)+"]", bbox={'facecolor': '0.85', 'pad': 10})

# Opens the plot
show()



