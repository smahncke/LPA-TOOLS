import numpy as np
from scipy.integrate import odeint
from pylab import *
import scipy.constants as const
from numba import double
from numba.decorators import jit, autojit
import scipy.constants as const
import scipy.integrate as integrate

def gaussian_field(max_a, z, r, t, w0, ctau, zf, lambda_0):


	k_0 = 2*np.pi/lambda_0

	# Calculate the Rayleigh length

	zr = 0.5*k_0*w0**2
	inv_zr = 1./zr
	inv_ctau2 = 1./ctau**2
    
    
	# Diffraction and stretch_factor
	diffract_factor = 1. - 1j*(z-zf)*inv_zr
	stretch_factor = 1.
	# Calculate the argument of the complex exponential
	exp_argument = 1j*k_0*( const.c*t + 0 - z )         - r**2 / (w0**2 * diffract_factor)         - 1./stretch_factor * inv_ctau2 * ( const.c*t  + 0 - z )**2
	# Get the transverse profile
	profile_Eperp = np.exp(exp_argument)         / ( diffract_factor * stretch_factor**0.5 )
    
	return(max_a*profile_Eperp.real )


def gaussian_envelope(max_a,z,ctau):
	rms = ctau/np.sqrt(2)
	return max_a*np.exp(-(z-0)**2/(2*rms**2))
