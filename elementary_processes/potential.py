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
from elementary_processes import ionization as ion
from particles import particle_data as ptcl
from particles import density_profiles as dens

#----------
# WAKEFIELD
#----------

#Plasma frequency
def wp(z,n_e): 
	"""
	Calculates the plasma frequency
	
	-INPUT:		- z [array]:	Array over the plasma target
			- n_e [float]:	Max. plasma density (in 1/cm^3)
	
	-RETURN: 	- wp [array]:	Plasma frequency at position z
	"""
	
	#Get the plasma density at position z
	n_gas = dens.homogen_density(z,n_e)
	
	#Return the plasma frequency
	return sqrt(n_gas*const.e**2/(const.epsilon_0*const.m_e))

#Plasma wavenumber
def kp(z,n_e):
	"""
	Calculates the plasma wavenumber out of the plasma frequency
	
	-INPUT:		- z [array]:	Array over the plasma target
			- n_e [float]:	Max. plasma density (in 1/cm^3)
			
	-RETURN:	- kp [array]:	Plasma wavenumber at position z
	"""
	
	return wp(z,n_e)/const.c

#Potential equation
def potentialeq(phi,zz,z,zi,max_a,ctau,n_e): 
	"""
	Non-linear, 1d differential wave-equation
	
	-INPUT:		-phi [array]:	Potential of the wake
			-zz [array]:	Array over the plasma target
			-zi [integer]:	Order of the solution of the wave-equation
			-max_a [float]:	Peak a0 of the laser
			-ctau [float]:	Laser duration (in femtoseconds)
			-n_e [float]:	Max. plasma density (in 1/cm^3)
	
	-RETURN:	-potentialeq [array]:	Wave-equation of the non-linear, 1D plasma wave
	"""
	
	return array([phi[1],((1+las.gaussian_envelope(max_a,zz,ctau)**2)/(2*(phi[0]+1))-0.5)*(kp(z,n_e)[zi]/1e6)**2])

#Potential/Solution of the wave-equation
def potential(zz,max_a,ctau,n_e):
	"""
	Gives the solution of the non-linear, 1d wave-equation
	
	-INPUT:		-zz [array]:	Array over the plasma target
			-max_a [float]:	Peak a0 of the laser
			-ctau [float]:	Laser duration (in femtoseconds)
			-n_e [float]:	Max. plasma density (in 1/cm^3)
	
	-RETURN:	-phi [array]:	Solution of the non-linear, 1d wave-equation
	"""
	#Order of the solution
	z = np.linspace(0, 1000, 32)
	
	#Initial conditions
	phiinit = array([0,0]) 
	
	#Solves the equation
	phi = list(map(lambda x: odeint(potentialeq, phiinit, zz, args = (z,x,max_a,ctau,n_e,)),range(10)))
	
	return phi
	
def wakefield(zz,max_a,ctau,n_e):
	"""
	Gives the wakefield 
	
	-INPUT:		-zz [array]:	Array over the plasma target
			-max_a [float]:	Peak a0 of the laser
			-ctau [float]:	Laser duration (in femtoseconds)
			-n_e [float]:	Max. plasma density (in 1/cm^3)
	
	-RETURN:	-wake [array]:	Wakefield
	"""
	#Order of the solution of the non-linear wave-equation
	z = np.linspace(0, 1000, 32)
	
	#Initial conditions
	phiinit = array([0,0]) 
	
	#Solves the non-linear wave-equation
	phi = list(map(lambda x: odeint(potentialeq, phiinit, zz, args = (z,x,max_a,ctau,n_e,)),range(10)))
    
	return(np.gradient(phi[3][:,0],70/1e6))


#-------------------
# TRAPPING CONDITION
#-------------------
	
def gammap(z,max_a,lambda_0): 
	"""
	Calculates gamma of the plasma
	
	-INPUT:		-z [array]:		Array over the plasma target
			-max_a [float]:		Peak a0 of the laser
			-lambda_0 [float]:	Laser wavelength (in microns)
	
	-RETURN:	-gammap [array]:	Gamma of the plasma
	"""
	#Frequency of the laser
	w_L = 2 * const.pi * const.c / (lambda_0*1e-6) 
	
	return w_L / wp(z) * sqrt((sqrt(1+max_a**2)+1)/2)
	
def get_phi_min(z,max_a,lambda_0):
	"""
	Calculates the minimum of the potential
	
	-INPUT:		-z [array]:		Array over the plasma target
			-max_a [float]:		Peak a0 of the laser
			-lambda_0 [float]:	Laser wavelength (in microns)
		
	-RETURN:	-phi_min [array]:	The minimum of the potential
	"""
	#The max of the potential 
	E_max = 2	
	
	return (-1+1/(2*gammap(z,max_a,lambda_0)**2))*(1-1/E_max**2) + E_max**2/(4*gammap(z,max_a,lambda_0)**2)


def condition(z,max_a,ctau,lambda_0,n_e):
	"""
	Calculates the condition for trapping (i.e. H_e < H_s, where H_e = Hamiltonian of the electron and 
	H_s = Hamiltonian of the separatrix)
	
	-INPUT:		-z [array]:		Array over the plasma target
			-max_a [float]:		Peak a0 of the laser
			-ctau [float]:		Laser duration (in femtoseconds)
			-lambda_0 [float]:	Laser wavelength (in microns)
			-n_e [float]:		Max. plasma density (in 1/cm^3)
			
	-RETURN:	cond [array]:		H_s - H_e
	"""
	#Calculate the solution of the wave-equation
	phi = potential(z,max_a,ctau,n_e)
	
	#Calculate the hamiltonian of the separatrix
	right_side = np.sqrt(1+las.gaussian_envelope(max_a,z,ctau)**2)/gammap(z,max_a,lambda_0)
	
	#Calculate the hamiltonian of the separatrix
	left_side =  1 + get_phi_min(z,max_a,lambda_0) - phi[3][:,0]
	
	return(right_side-left_side)
	
def condition_fullfilled(z,max_a,ctau,lambda_0,n_e):
	"""
	Gives the z positions of the target where trapping can happen
	
	-INPUT:		-z [array]:		Array over the target
			-max_a [float]:		Peak a0 of the laser
			-ctau [float]:		Laser duration (in femtoseconds)
			-lambda_0 [float]:	Laser wavelength (in microns)
			-n_e [float]:		Max. plasma density (in 1/cm^3)
	
	-RETURN:	cond_fullfill. [array]:	Positions at the target, where H_e < H_s (trapping condition fullfilled)
	"""
	#Calculate the trapping condition first
	cond = condition(z,max_a,ctau,lambda_0,n_e)
	
	#Initialize the array
	fullfilled = np.zeros_like(z)
	
	#If H_s > H_e the ionized electrons get trapped, otherwise they don't
	for i in range(0,len(z)):
		if cond[i] > 0:
			fullfilled[i] = cond[i]
		else:
			fullfilled[i] = 0
	return(fullfilled)
