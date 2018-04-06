import numpy as np
from scipy.integrate import odeint
from pylab import *
import scipy.constants as const
from numba import double
from numba.decorators import jit, autojit
import scipy.constants as const
import scipy.integrate as integrate

#-------------------------------
# Ionization energy of the atoms
#-------------------------------

def get_ion_energy(element,niveau):

	"""
	Returns the ionization of the energy of the given element and niveau
	
	-INPUT:	- element [string]:	Name of the element (shortcut, i.e. 'N' for nitrogen)
		- niveau [integer]:	Ionization niveau/level
		
	-RETURN: - U_i [float]:		Ionization energy of the element and level (in eV)
	"""
	
	#Get the energy and return an error if the element doesn't exist
	if element == 'N':
		U_G = [14.53,29.60,47.44,77.47,97.89,552.07,667.05]
	elif element == 'Ar':
		U_G = [15.76,27.62,40.74,59.81,75.02,91.01,124.32,\
				143.46,422.45,478.69,538.96,618.26,686.10,\
				755.74,854.77,918.03,4120.885,4426.23]
	elif element == 'Kr':
		U_G = [14.0,24.36,36.95,52.5,64.7,78.5,111.0,125.80,\
				230.85,268.2,308.,350.,391.,447.,492.,541.,\
				592.,641.,786.,833.,884.,937.,998.,1051.,1151.,\
				1205.3,2928.,3070.,3227.,3381.]
	else:
		print("Error: Unknown element given!")
	
	#Returns the energy for total ionization if niveau > number of electrons of the atom
	if niveau >= len(U_G):
		U_i = U_G[len(U_G)-1]
	else:
		U_i = U_G[niveau-1]

	return(U_i)
