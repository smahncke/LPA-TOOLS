#IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from pylab import title, show
from scipy.constants import e, c, m_e, h, alpha, hbar

#Initial parameters
lambda_0 = 0.8e-6 #laser wavelength
a = 1.5 #Max laser amplitude
E_L = 8.e12# Field of the laser

#---------------------------------------
# Don't change anything below this line
#---------------------------------------

#Physical constants
lambda_c = h/(m_e * c) #compton wavelength
U_H = 13.6 #ionization potential of hydrogen in eV

#Other parameters
k_0 = 1/lambda_0 #wave number of the laser
E_k = (1/e)*k_0*(c**2)*m_e #laser energy
omega_0 = c*k_0 #laser frequency

#ionization probability
def ion_prob(epsilon,U_i):
    """ Inputs: 
        - epsilon (1D array): Energy of the ionized electron
        - U_i (float): Ionization potential 
        Output: Probability, that the electron got ionized with an 
        initial energy epsilon
        """
    #Keldysh parameter
    gamma_k = (alpha/a)*np.sqrt(U_i/U_H) 
    
    #Argument of the exponential function
    laser_ion = lambda_0/lambda_c * a**3 *gamma_k**3 * (E_k/(E_L))
    tunnel_ion = (gamma_k**3 * epsilon)/(hbar*omega_0)
    
    #Probability
    prob = np.exp(-2/3 * (laser_ion + tunnel_ion))
    
    return(prob)

#The possible energy range
energy = np.linspace(0,e*1.e6,1000)

### Plot the functions

# Set the size of the plot
fig, ax = plt.subplots(figsize=(10,10))

# Set the labels of the axis
plt.xlabel("Initial energy $\epsilon$ [eV]")
plt.ylabel(r"Probability $w(\epsilon)$")

ax.grid(which='major', alpha=0.8, color='grey', linestyle='--', linewidth=0.6, animated='True')
ax.grid(which='minor', alpha=0.6, color='grey', linestyle='--', linewidth=0.5, animated='True')

# Make the plots
ax.plot(energy, ion_prob(energy,77.472), color='grey', label="$N^{4+}$")
ax.plot(energy, ion_prob(energy,97.888), color='blue', label="$N^{5+}$")
ax.plot(energy, ion_prob(energy,552.057), color='green', label="$N^{6+}$")
ax.plot(energy, ion_prob(energy,667.029), color='red', label="$N^{7+}$")

# Make the legend 
legend = ax.legend(loc='upper right', shadow=True)
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
title("Ionization probability [ $\lambda_0 =$ "+str(lambda_0*1.e9)+" nm, $a_0 =$ "+str(a)+", $|E_L| =$"+str(E_L/1.e12)+"TV/m]", bbox={'facecolor': '0.85', 'pad': 10})

# Opens the plot
show()


