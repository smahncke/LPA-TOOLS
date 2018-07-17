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

def gaussian_field_OLD(z, max_a, w0, ctau, zf, lambda_0, r = 0, t = 0):

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

#-------------------
# Gaussian Field NEW
#-------------------

def gaussian_field(z,a0,lambda_0,w0,ctau,zf,x=0, y=0, t=0,theta_pol=0, phi2_chirp=0,cep_phase=0,z0=0):
    """
        See the docstring of LaserProfile.E_field
        """
    # Note: this formula is expressed with complex numbers for compactness
    # and simplicity, but only the real part is used in the end
    # (see final return statement)
    # The formula for the laser (in complex numbers) is obtained by
    # multiplying the Fourier transform of the laser at focus
    # E(k_x,k_y,\omega) = exp( -(\omega-\omega_0)^2(\tau^2/4 + \phi^(2)/2)
    # - (k_x^2 + k_y^2)w_0^2/4 ) by the paraxial propagator
    # e^(-i(\omega/c - (k_x^2 +k_y^2)/2k0)(z-z_foc))
    # and then by taking the inverse Fourier transform in x, y, and t
    
    # Diffraction and stretch_factor
    
    k0 = 2*np.pi/lambda_0
    E0 = a0*const.m_e*const.c**2*k0/const.e
    zr = 0.5*k0*w0**2

    inv_zr = 1./zr
    E0x = E0 * np.cos(theta_pol)
    E0y = E0 * np.sin(theta_pol)
    inv_ctau2 = 1./(ctau)**2
    
    prop_dir = 1
    diffract_factor = 1. + 1j * prop_dir*(z - zf) * inv_zr
    stretch_factor = 1 - 2j * phi2_chirp * const.c**2 * inv_ctau2
    
    # Calculate the argument of the complex exponential
    exp_argument = - 1j*cep_phase \
        + 1j*k0*( prop_dir*(z - z0) - const.c*t ) \
        - (x**2 + y**2) / (w0**2 * diffract_factor) \
        - 1./stretch_factor*inv_ctau2 * \
            ( prop_dir*(z - z0) - const.c*t )**2
    # Get the transverse profile
    profile = np.exp(exp_argument) /(diffract_factor * stretch_factor**0.5)
    
    # Get the projection along x and y, with the correct polarization
    Ex = E0x * profile
    Ey = E0y * profile
    
    return( Ex.real)
