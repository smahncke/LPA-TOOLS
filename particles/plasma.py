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

def get_critical_density(lambda_0):
    
    epsilon_0 = 8.854e-12

    omega_0 = 2*np.pi*const.c/lambda_0

    n_c = omega_0**2 * epsilon_0 * const.m_e/const.e**2

    return(n_c)
