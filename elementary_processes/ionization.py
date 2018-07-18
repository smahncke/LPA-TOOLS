import numpy as np
from scipy.integrate import odeint
from pylab import *
import scipy.constants as const
from numba import double
from numba.decorators import jit, autojit
import warnings
import scipy.constants as const
import scipy.integrate as integrate
import sys
sys.path.append('../')
from laser import laser_profiles as las
from sklearn.preprocessing import normalize

#-----------------------
# Ionization probability
#-----------------------


def ionization_probability(z,U_i,max_a,ctau,lambda_0=0.8,energy=0):

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
    E_L = las.gaussian_field_envelope(z, max_a,ctau,lambda_0) 
    
	#Envelope of the laser
    amplitude = las.gaussian_a0(z, max_a,ctau)


	#Keldysh parameter
    gamma_k = (const.alpha/amplitude)*np.sqrt(U_i/U_H) 
    
    #Argument of the exponential function
    laser_ion = lambda_0/lambda_c * amplitude**3 *gamma_k**3 * (E_k/(np.abs(E_L)))
    tunnel_ion = (gamma_k**3 * energy)/(const.hbar*omega_0)
    
   	#ADK rate
    rate = np.exp(-2/3 * (laser_ion + tunnel_ion))

    return(rate)

#---------------------
# Degree of ionization
#---------------------

def ionization_degree(z,probability):
	
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
		
	return(degree/ne)


def get_min_a0(z,U_i,ctau,max_a = 10,ion_thres = 0.01):

	"""
	Calculates the lowest value of a0 that is needed for ionization

	-INPUT:	-z [array]:				Array over the laser pulse
			-w0 [float]:			Waist of the laser in the focus (in microns)
			-ctau [float]:			Laser duration (in femtoseconds)
			-zf [float]:			Laser focus (in microns)
			-lambda_0 [float]:		Laser wavelength (in microns)
			-U_i [float]:			Ionization energy (in eV)
			-energy [float]:		Final energy of the ionized electrons (in J)
			-max_a [float,opt]:		The upper limit for the laser peak a0 
			-ion_thres [float,opt]:	The threshold/ minimum ionization degree

	-RETURN: -min_a0 [float]:		The lower limit of a0 for ionization
			
	"""
	#Calculate the ionization probability
	prob = ionization_probability(z,U_i,max_a,ctau)

	#Calculate the ionization degree
	degree = ionization_degree(z,prob)

	#Get the laser envelope
	envelope = las.gaussian_a0(z,max_a,ctau)

	#Check for which value of a0 the ionization degree is higher than the threshold
	for i in range(0,len(degree)):
		if degree[i] > ion_thres:
			min_a0 = envelope[i]
			break
	
	return(min_a0)	

