import numpy as np
from scipy.integrate import odeint
from pylab import *
import scipy.constants as const
from numba import double
from numba.decorators import jit, autojit
import scipy.constants as const
import scipy.integrate as integrate
import sys
sys.path.append('../')
from laser import laser_profiles as las

#-----------------------
# Ionization probability
#-----------------------

def ionization_probability(z,max_a,w0,ctau,zf,lambda_0,U_i,energy, plot=False, r = 0, t = 0):

	"""
	Calculates the ionization probability in dependence of laser parameters
	
	INPUT:	- U_i [float]: 		Ionization energy of the atom and the ionization level
		- energy [float]:	Kinetic energy of the ionized electron after ionization (in MeV)
		- max_a [float]:	Peak a0 of the laser
		- z [array]:		Array over the laser pulse
		- w0 [float]:		Waist of the laser in the focus (in microns)	
		- ctau [float]:		Laser duration (in femtoseconds)
		- zf [float]:		Focus position of the laser (in microns)
		- lambda_0 [float]:	Laser wavelength (in microns)
		- plot [boolean,opt.]:	Wether to plot the result
		- r [float, opt.]:	Radial position of the laser (in microns)
		- t [float, opt.]:	Temporal position of the laser
		
	RETURN:	- prob [array]: Ionization probability
	"""
	#Constants and important parameters
	U_H = 13.6
	epsilon_0 = 8.854e-12
	lambda_c = const.h/(const.m_e * const.c)
	omega_0 = 2*np.pi*const.c/lambda_0
	E_k = (1/const.e)*(2*np.pi/lambda_0)*(const.c**2)*const.m_e
    
	#Electrical field of the laser
	E_gauss = las.gaussian_field(z, max_a,w0,ctau,zf,lambda_0) 
    
	#Envelope of the laser
	amplitude = las.gaussian_envelope(z, max_a,ctau)

	#Keldysh parameter
	gamma_k = (const.alpha/amplitude)*np.sqrt(U_i/U_H) 
    
	#Calculate the electrical field of the laser (in V/m)
	E_L = 1.8e2*np.sqrt(2.8e18/(lambda_0*const.c*epsilon_0))*E_gauss #Electrical field of the laser in V/m
	
    	#Argument of the exponential function
	laser_ion = lambda_0/lambda_c * amplitude**3 *gamma_k**3 * (E_k/(np.sqrt(E_L**2)))
	tunnel_ion = (gamma_k**3 * energy)/(const.hbar*omega_0)
    
   	#Probability
	prob = np.exp(-2/3 * (laser_ion + tunnel_ion))

	#Returns the plot if 'plot' is True
	if plot == True:
		plt.figure()
		plt.plot(z,E_gauss)
		plt.plot(z,prob)
		plt.show()

	return(prob)

#---------------------
# Degree of ionization
#---------------------

def ionization_degree(z,probability,plot = False):
	
	"""
	Calculates the degree of ionization behind the laser pulse
	
	INPUT:	- probability [array]:	Ionization probability of the gas in dependence of the laser and gas parameters
		- z [array]:		Array over the laser pulse, should have the same dimension as the probability
		- plot [boolean, opt.]:	Wether to plot the result
		
	RETURN:	- degree/ne [array]:	Degree of ionization at the position z in percent
	"""
	
	#Normalize the ionization degree
	ne = 100
	
	#Initialize the array
	degree = np.zeros_like(z)
    
	#Calculate the degree of ionization
	for i in range(0,len(z)):
		if i+1 < len(z):
			degree[i+1] = (degree[i]+probability[i]*(ne-degree[i]))
	
	#Plot the result of 'plot' is True
	if plot == True:
		plt.figure()
		plt.plot(z,probability)
		plt.plot(z,degree/ne)
		plt.show()
		
	return(degree/ne)
