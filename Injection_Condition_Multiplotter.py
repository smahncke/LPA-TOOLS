
# coding: utf-8

# In[1]:

import numpy as np
from scipy.integrate import odeint
from pylab import *
import scipy.constants as const
from numba import double
from numba.decorators import jit, autojit
from scipy.constants import m_e,c,e, epsilon_0
import matplotlib.pyplot as plt
from pylab import title, show



# #### Input

# In[250]:

n_e = 1.e24 # Plasma density
a0 = 0. # normalized laser amplitude
lambda_0 = 800.e-9 #laser wavelength
gamma_p = 33
RMS = 1
E_max = 2



# ##### Parameters

# In[251]:

omega_0 = 2*np.pi*c/lambda_0


# ##### Potential

# In[252]:

def n(zz):  # Density distribution
    # return where(zz>500, ne/2. * ones(zz.size), ne * np.ones(size(zz)))
    return n_e * ones(zz.size)

# Plasma frequency
def omega_p(zz): 
    return sqrt(n(zz)*e**2/(epsilon_0*m_e))

 # Plasma phase velocity
def beta_p(zz):
    return 1-omega_p(zz)**2/(2*omega_0**2)

# Plasma wavenumber
def k_p(zz): 
    return omega_p(zz)/c

# Laser profile
def a(psi,a0,RMS):
    return a0*np.exp(-psi**2/(4*RMS**2))

# Potential ORIGINAL
#def potentialeq(phi,xi,z,zi): # Potential equation
#    return array([phi[1],((1+a(xi*k_p(z)[zi],a0,RMS)**2)/(2*(phi[0]+1))-0.5)*(kp(z)[zi]/1e6)**2])

# Potential NEW
def potentialeq(phi,psi,z,zi):
    return array([phi[1],gamma_p**2 * (beta_p(z)[zi]*(np.sqrt(1/(1 - (1+ a(psi,a0,RMS)**2)/(gamma_p**2 * (1+phi[0])**2))))-1)])

def get_phi_min(E_max):
    return (-1+1/(2*gamma_p**2))*(1-1/E_max**2) + E_max**2/(4*gamma_p**2)


# ##### Plots

# In[253]:



fig, ax = plt.subplots(figsize=(12,5))
#ax.plot(psi/np.pi,left_side())
#ax.plot(psi/np.pi,right_side())
for i in range(0,10):
	psi = np.linspace(-2*np.pi,2*np.pi,1e2)
	z = np.linspace(0, 1000, 100)

#xi = linspace(0,-200.0,100)	# Range to solve phi

	phiinit = array([0,0]) #initial conditions
	phi = list(map(lambda x: odeint(potentialeq, phiinit, psi, args = (z,x,)),range(10)))

#fig, ax = plt.subplots(figsize=(12,5))
#ax.plot(psi/np.pi,phi[3][:,0])
#ax.plot(xi,nwake(phi,xi,3)-1)
#ax.plot(psi/np.pi,10*a(psi,a0,RMS))
#ax.plot(xi,wake(phi,xi,3)*5)
#plot(xi,5*sin(xi*kp(z)[3]/1e6))
#ax.set_xlabel(r'$\psi$',fontdict={'fontsize':20})
#ax.set_ylabel('$\phi$',fontdict={'fontsize':20})
#ax.set_ylim(-10,20);


# ##### Trapping condition

# In[ ]:
#	def get_phi_min(a0,RMS):
#		psi_min = int(25*(np.sqrt(4*RMS**2 * np.log(a0/1e-3))/np.pi))		
#		phi_min = phi[3][psi_min,0]
#		return(phi_min)


# LEFT SIDE
	def left_side():
		return 1 + get_phi_min(E_max) - phi[3][:,0]

# Right SIDE
	def right_side():
		return np.sqrt(1+a(psi,a0,RMS)**2)/gamma_p

# Condition
	def condition_fullfilled():
		diff = np.zeros_like(psi, dtype=object)
		diff = np.where(right_side()-left_side() < 0,0,right_side()-left_side())
		return diff

	ax.plot(psi/np.pi,condition_fullfilled(),label="$a_0 =$ "+str("%.1f" % a0)+"")
	#ax.plot(psi/np.pi,a(psi,a0,RMS), linestyle='--')
	a0 += 0.2

ax.plot(psi/np.pi,a(psi,1.7,RMS),linestyle='--', color="grey", label="Laser")
legend = ax.legend(loc='upper left', shadow=True)
for label in legend.get_texts():
    label.set_fontsize('medium')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

ax.set_xlabel(r'$\psi\,[\pi]$',fontdict={'fontsize':20})
ax.set_ylabel('$H_S - H_i$',fontdict={'fontsize':20})
show()


# In[ ]:



