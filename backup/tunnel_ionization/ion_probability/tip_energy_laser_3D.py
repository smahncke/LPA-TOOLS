"""
- RESPOSITORY:	LPA-II TOOLS
- CATEGORY: 	TUNNELING IONIZATION PROBABILITY
- DATE:			10/30/2017

- INFO:	Plots a 3D plot that shows the probability of tunneling ionization in dependence of
		the initial electron energy after ionization and the laser amplitude a for 
		a choosable element.		

"""

#IMPORTS
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from pylab import title, show
from scipy.constants import e, c, m_e, h, alpha, hbar

#Initial parameters
lambda_0 = 800e-9 #laser wavelength in m
element = 'Ar' # Element that is used fo ionization injection (string, use the shortcut, f.e. 'N' for nitrogen)
ionization_level = 5 #The atom level that will be ionized (f.e. 5 means: 5+ -> 6+)


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

#ionization probability
def ion_prob(epsilon,U_i,amplitude):
	""" Inputs: 
        - epsilon (1D array): Energy of the ionized electron
        - U_i (float): Ionization potential 
		- amplitude (1D array): laser amplitude
        Output: Probability, that the electron got ionized with an 
        initial energy epsilon (float)
	"""
	#Keldysh parameter
	gamma_k = (alpha/amplitude)*np.sqrt(U_i/U_H) 
    
	#Calculate the electric field of the laser in V/m
	E_L = 1e2*np.sqrt(2.8e18/(lambda_0*1e6*c*epsilon_0))*amplitude 

    #Argument of the exponential function
	laser_ion = lambda_0/lambda_c * amplitude**3 *gamma_k**3 * (E_k/(E_L))
	tunnel_ion = (gamma_k**3 * epsilon)/(hbar*omega_0)
    
    #Probability
	prob = np.exp(-2/3 * (laser_ion + tunnel_ion))
    
	return(prob)

#The possible energy and laser amplitude range
energy = np.linspace(0,e*1e6,2500)
a0 = np.linspace(0.01,2,2500)
energy, a0 = np.meshgrid(energy, a0)

### Plot the functions

fig = plt.figure()
ax = fig.gca(projection='3d')

# Set the labels of the axis
plt.xlabel("Initial electron energy $\epsilon$ [MeV]")
plt.ylabel("Laser amplitude a")
ax.set_zlabel(r"Probability $w(\epsilon,a)$")


#Create the plot
surf = ax.plot_surface(energy/(e*1e6), a0, ion_prob(energy,U_G[ionization_level],a0), cmap=cm.plasma,
                       linewidth=0, antialiased=False)

#Customize the z-axis
ax.set_zlim(0,1.00)

# Make the legend 

fig.colorbar(surf, shrink=0.5, aspect=5)

# Open the plot
show()



