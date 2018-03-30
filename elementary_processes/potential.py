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

def wp(z,n_e): # Plasma frequency
	n_gas = dens.homogen_density(z,n_e)
	return sqrt(n_gas*const.e**2/(const.epsilon_0*const.m_e))

def kp(z,n_e): # Plasma wavenumber
    return wp(z,n_e)/const.c
	
def potentialeq(phi,zz,z,zi,max_a,ctau,n_e): # Potential equation
    return array([phi[1],((1+las.gaussian_envelope(max_a,zz,ctau)**2)/(2*(phi[0]+1))-0.5)*(kp(z,n_e)[zi]/1e6)**2])
	
def potential(zz,max_a,ctau,n_e):
	z = np.linspace(0, 1000, 32)
	phiinit = array([0,0]) #initial conditions
	phi = list(map(lambda x: odeint(potentialeq, phiinit, zz, args = (z,x,max_a,ctau,n_e,)),range(10)))
	return phi
	
def wakefield(zz,max_a,ctau,n_e):
	z = np.linspace(0, 1000, 32)
	phiinit = array([0,0]) #initial conditions
	phi = list(map(lambda x: odeint(potentialeq, phiinit, zz, args = (z,x,max_a,ctau,n_e,)),range(10)))
    
	return(np.gradient(phi[3][:,0],70/1e6))

#-------------------
# TRAPPING CONDITION
#-------------------
	
def gammap(z,max_a,w0,lambda_0): 
    w_L = 2 * const.pi * const.c / (lambda_0*1e-6) 
    return w_L / wp(z) * sqrt((sqrt(1+max_a**2)+1)/2)
	
def get_phi_min(z,max_a,ctau,lambda_0):
	E_max = 2
    return (-1+1/(2*gammap(z,max_a,ctau,lambda_0)**2))*(1-1/E_max**2) + E_max**2/(4*gammap(z,max_a,ctau,lambda_0)**2)


def condition(z,max_a,ctau,lambda_0,n_e):
	phi = potential(z,max_a,ctau,n_e)
    right_side = np.sqrt(1+las.gaussian_envelope(max_a,z,ctau)**2)/gammap(z,max_a,w0,lambda_0)
    left_side =  1 + get_phi_min(z,max_a,ctau,lambda_0) - phi[3][:,0]
    return(right_side-left_side)
	
def condition_fullfilled(z,max_a,ctau,lambda_0,n_e):
    cond = condition(z,max_a,ctau,lambda_0,n_e)
    fullfilled = np.zeros_like(z)
    for i in range(0,len(z)):
        if cond[i] > 0:
            fullfilled[i] = cond[i]
        else:
            fullfilled[i] = 0
    return(fullfilled)