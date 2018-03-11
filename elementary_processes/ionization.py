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



def ionization_probability(U_i,energy,max_a,z,r,t,w0,ctau,zf,lambda_0, plot=False):

	U_H = 13.6
	epsilon_0 = 8.854e-12
	lambda_c = const.h/(const.m_e * const.c)
	omega_0 = 2*np.pi*const.c/lambda_0
	E_k = (1/const.e)*(2*np.pi/lambda_0)*(const.c**2)*const.m_e
    
	E_gauss = las.gaussian_field(max_a, z,r,t,w0,ctau,zf,lambda_0) 
    
	amplitude = las.gaussian_envelope(max_a,z,ctau)

	#Keldysh parameter
	gamma_k = (const.alpha/amplitude)*np.sqrt(U_i/U_H) 
    
	#Calculate the electrical field of the laser (in V/m)
	E_L = 1.8e2*np.sqrt(2.8e18/(lambda_0*const.c*epsilon_0))*E_gauss #Electrical field of the laser in V/m
	
    #Argument of the exponential function
	laser_ion = lambda_0/lambda_c * amplitude**3 *gamma_k**3 * (E_k/(np.sqrt(E_L**2)))
	tunnel_ion = (gamma_k**3 * energy)/(const.hbar*omega_0)
    
    #Probability
	prob = np.exp(-2/3 * (laser_ion + tunnel_ion))

	if plot == True:
		plt.figure()
		plt.plot(z,E_gauss)
		plt.plot(z,prob)
		plt.show()

	return(prob)

def ionization_degree(probability,z,plot = False):
	
	ne = 100
	degree = np.zeros_like(z)
    

	for i in range(0,len(z)):
		if i+1 < len(z):
			degree[i+1] = (degree[i]+probability[i]*(ne-degree[i]))

	if plot == True:
		plt.figure()
		plt.plot(z,probability)
		plt.plot(z,degree/ne)
		plt.show()
	return(degree[len(z)-1]/ne,degree/ne)
