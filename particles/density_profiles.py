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

def homogen_density(z,max_n):
	return max_n*np.ones_like(z)
