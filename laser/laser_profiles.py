import numpy as np
from scipy.integrate import odeint
from pylab import *
import scipy.constants as const
from numba import double
from numba.decorators import jit, autojit
import scipy.constants as const
import scipy.integrate as integrate

#------------------------------
# Electrical field of the laser
#------------------------------

def gaussian_field(z, max_a, w0, ctau, zf, lambda_0, r = 0, t = 0):

	"""
	Calculates the electrical field of the gaussian laser pulse
	
	-INPUT:	- max_a [float]:	Peak a0 of the laser pulse
		- z [array]:		Array over the laser pulse
		- w0 [float]:		Waist of the laser pulse in the focus (in microns)
		- ctau [float]:		Laser duration (in femtoseconds)
		- zf [float]:		Laser position (in microns)
		- lambda_0 [float]:	Laser wavelength (in microns)
		- r [float]:		Radial position of the laser pulse (in microns)
		- t [float]:		Temporal position of the laser pulse
		
	-RETURN: - field [array]:	Electrical field of the laser pulse
	"""
	
	#Wave number
	k_0 = 2*np.pi/lambda_0

	# Calculate the Rayleigh length
	zr = 0.5*k_0*w0**2
	inv_zr = 1./zr
	inv_ctau2 = 1./ctau**2
    
    
	# Diffraction and stretch_factor
	diffract_factor = 1. - 1j*(z-zf)*inv_zr
	stretch_factor = 1.
	
	# Calculate the argument of the complex exponential
	exp_argument = 1j*k_0*( const.c*t + 0 - z ) - r**2 / (w0**2 * diffract_factor) - 1./stretch_factor * inv_ctau2 * ( const.c*t  + 0 - z )**2
	
	# Get the transverse profile
	profile_Eperp = np.exp(exp_argument)         / ( diffract_factor * stretch_factor**0.5 )
    
	return(max_a*profile_Eperp.real )

#-------------------------------------
# Envelope of the gaussian laser pulse
#-------------------------------------

def gaussian_envelope(z, max_a, ctau):
	
	"""
	Calculates the envelope of the gaussian laser pulse
	
	-INPUT: - max_a [float]:	Peak a0 of the laser pulse
		- z [array]:		Array over the laser pulse
		- ctau [float]:		Laser duration (in femtoseconds)
		
	-RETURN: -envelope [array]:	Envelope of the gaussian laser pulse
	"""
	#Calculates the root-mean-square of the laser pulse
	rms = ctau/np.sqrt(2)
	
	return max_a*np.exp(-(z-0)**2/(2*rms**2))
