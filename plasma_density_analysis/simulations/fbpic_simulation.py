"""
This is an advanced FBPIC input script that uses the PDA script for
gas density profile calculation.

Usage
-----
- Modify the parameters below to suit your needs
- Type "python -i lwfa_script.py" in a terminal
- When the simulation finishes, the python session will *not* quit.
    Therefore the simulation can be continued by running sim.step()

Help
----
All the structures implemented in FBPIC are internally documented.
Enter "print(fbpic_object.__doc__)" to have access to this documentation,
where fbpic_object is any of the objects or function of FBPIC

Visit the https://github.com/fbpic/fbpic or https://github.com/smahncke/PDA 
for further information.
"""

# -------
# Imports
# -------
import numpy as np
from scipy.constants import c, m_p, m_e, e
# Import the relevant structures in FBPIC
from fbpic.main import Simulation, adapt_to_grid
from fbpic.particles import Particles
from fbpic.lpa_utils.laser import add_laser
from fbpic.openpmd_diag import FieldDiagnostic, ParticleDiagnostic, \
     set_periodic_checkpoint, restart_from_checkpoint


# ----------
# Parameters
# ----------

# Whether to use the GPU
use_cuda = True

# The simulation box
Nz = 4500         # Number of gridpoints along z 
zmax = 30.e-6    # Right end of the simulation box (meters)
zmin = -60.e-6   # Left end of the simulation box (meters)
Nr = 300          # Number of gridpoints along r
rmax = 100.e-6    # Length of the box along r (meters)
Nm = 2           # Number of modes used

# The simulation timestep
dt = (zmax-zmin)/Nz/c   # Timestep (seconds)
N_step = int(5.1e-3/(dt*c))    # Number of iterations to perform #ORIGINAL

# Order of accuracy of the spectral, Maxwell (PSATD) solver.
# -1 correspond to infinite order, i.e. wave propagation is perfectly
# dispersion-free in all directions. This is adviced for single GPU/CPU
# simulations. For multi GPU/CPU simulations, choose n_order > 4
# (and multiple of 2). A large n_order leads to more overhead in MPI
# communications, but also to a more accurate dispersion for waves.
# (Typically, n_order = 32 gives accurate physical results)
n_order = -1

# The particles
p_zmin = 25.e-6  # Position of the beginning of the plasma (meters)
p_zmax = 31.e-6  # Position of the end of the plasma (meters)
p_rmin = 0.      # Minimal radial position of the plasma (meters)
p_rmax = 90.e-6  # Maximal radial position of the plasma (meters)
p_nz = 2         # Number of electrons per cell along z
p_nr = 2         # Number of electrons per cell along r
p_nt = 4         # Number of electrons per cell along theta
N_nz = 2         # Number of nitrogen atoms per cell along z
N_nr = 2         # Number of nitrogen atoms per cell along r
N_nt = 4         # Number of nitrogen atoms per cell along theta


# The laser
a0 = 1.5          # Laser amplitude
w0 = 20.e-6       # Laser waist
ctau = 10.e-6     # Laser duration
z0 = 0.e-6      # Laser centroid
zf = 500.e-6    # Laser focus

# The moving window
v_window = c       # Speed of the window

# The diagnostics and the checkpoints/restarts
diag_period = int(N_step/100)         # Period of the diagnostics in number of timesteps
save_checkpoints = False # Whether to write checkpoint files
checkpoint_period = int(N_step/21)   # Period for writing the checkpoints
use_restart = False      # Whether to restart from a previous checkpoint
track_electrons = True  # Whether to track and write particle ids

#Gas and plasma properties
n_tot = 1.e18*1e6    #Total plasma density
am_N2 = 100.          #Amount of nitrogen (in percent)
FWHM = 300.e-6       #Full-width at half-max of the gaussian density profile (in m)

#Target properties 
ramp_start = 30.e-6      #Start point of the resulting plasma upramp (in m)
ramp_length = 500.e-6    #The length of the upramp and the downramp (in m)
plateau = 4000.e-6       #The length of the plateau (in m)

# Other settings
inj_thres = 6           #The ionization level of N2, at which the electron injection starts


#----------------
# DENSITIES
#----------------


# Gaussian profile of the nitrogen
mu = ramp_start + ramp_length           #mu of the gaussian nitrogen density profile
sigma = FWHM/(2*(np.sqrt(2*np.log(2))))  #sigma of the gaussian density profile

am_H2 = 100. - am_N2 #Amount of hydrogen

#The final plasma density function
def dens_func_plasma(z,r):
    
    #Initialize the array
    e_dens_tot = np.zeros_like(z, dtype=object)
    
    #Draw the plasma density under specific conditions
    e_dens_tot = np.where(z<ramp_start,0,e_dens_tot)
    e_dens_tot = np.where(np.logical_and(z<ramp_start+ramp_length, z>=ramp_start),n_tot*np.sin((0.5*np.pi*(z-ramp_start))/(ramp_length))**2, e_dens_tot)
    e_dens_tot = np.where(np.logical_and(z >= ramp_start+ramp_length, z <= ramp_start+ramp_length+plateau),n_tot,e_dens_tot)

    e_dens_tot = np.where(np.logical_and(z < ramp_start + 2*ramp_length + plateau,z > ramp_start+ramp_length+plateau),n_tot*np.sin((0.5*np.pi*(z-ramp_start))/ramp_length)**2, e_dens_tot)
    e_dens_tot = np.where(z>=ramp_start + 2*ramp_length + plateau,0,e_dens_tot)
    
    return(e_dens_tot)

#The gas density function of nitrogen 
def dens_func_N2(z,r):
    
    #Load the plasma density
    e_dens_tot = dens_func_plasma(z,r)
    
    #Initialize the arrays
    gas_dens_N2 = np.zeros_like(z,dtype=object)
    temp_dens = gas_dens_N2.copy()
    
    #Calculate the density
    temp_dens= 0.5*e_dens_tot*(1/(1-(am_N2/100 -(inj_thres-1)*(am_N2/100))*np.exp(-(z-mu)**2/(2*sigma**2))))
    gas_dens_N2 = (temp_dens-0.5*e_dens_tot)*(1/(2-inj_thres))*np.exp(-(z-mu)**2/(2*sigma**2))
    
    return(gas_dens_N2)

#The gas density function of hydrogen
def dens_func_H2(z,r):
    
    #Load the plasma and the nitrogen density
    e_dens_tot = dens_func_plasma(z,r)
    gas_dens_N2 = dens_func_N2(z,r)
    
    #Initialize the arrays
    gas_dens_H2 = np.zeros_like(z, dtype=object)
    e_dens_N2 = gas_dens_H2.copy()
    e_dens_H2 = gas_dens_H2.copy()
    
    # Electron density of the nitrogen
    e_dens_N2 = 2.*gas_dens_N2*(inj_thres-1)
    
    # Electron density of the hydrogen
    e_dens_H2 = e_dens_tot - e_dens_N2
    
    # Gas density of the hydrogen
    gas_dens_H2 = e_dens_H2/2
    
    return(gas_dens_H2)

# The electron density that causes the acceleration
def dens_func_e(z,r):
    
    # Load the gas densities
    gas_dens_N2 = dens_func_N2(z,r)
    gas_dens_H2 = dens_func_H2(z,r)

    #Calculate the electron density out of the gas densities
    electron_density = 2.*(inj_thres-1)*gas_dens_N2 + 2.*gas_dens_H2

    return(electron_density)

# The N5+ density
def dens_func_N2_thres(z,r):
    
    # Load the N2 gas density
    gas_dens_N2 = dens_func_N2(z,r)

    # Calculate the threshold out of the gas density
    N2_thres_dens = 2.*gas_dens_N2

    return(N2_thres_dens)

    

# ---------------------------
# Carrying out the simulation
# ---------------------------

# NB: The code below is only executed when running the script,
# (`python -i lpa_sim.py`), but not when importing it (`import lpa_sim`).

"""NOTE:    The density parameters n_H2 und n_N2 in the header of the 
            simulation particle functions should be 1, because the
            density functions already include these parameters"""

if __name__ == '__main__':


    # Initialize the simulation object
    sim = Simulation( Nz, zmax, Nr, rmax, Nm, dt,
        p_zmin, p_zmax, p_rmin, p_rmax, p_nz, p_nr, p_nt, 1,
        dens_func=dens_func_e, zmin=zmin, boundaries='open',
        n_order=n_order, use_cuda=use_cuda )
    
    p_zmin, p_zmax, Npz = adapt_to_grid( sim.fld.interp[0].z,
                                        p_zmin, p_zmax, N_nz )
    p_rmin, p_rmax, Npr = adapt_to_grid( sim.fld.interp[0].r,
                                        p_rmin, p_rmax, N_nr )
          
    # Add the non ionized nitrogen atoms 
    sim.ptcl.append(
        Particles(e, 14*m_p, 1, Npz, p_zmin,
            p_zmax, Npr, p_rmin, p_rmax,
            N_nt, dt, use_cuda=use_cuda, dens_func=dens_func_N2_thres,
            grid_shape=sim.fld.interp[0].Ez.shape,
            continuous_injection=True ) )
            
    sim.ptcl[1].make_ionizable(element='N', level_start=(inj_thres-1),
                               target_species=sim.ptcl[0])
                               
    # Load initial fields
    # Add a laser to the fields of the simulation
    add_laser( sim, a0, w0, ctau, z0, zf=zf )

    if use_restart is False:
        # Track electrons if required (species 0 correspond to the electrons)
        if track_electrons:
            sim.ptcl[0].track( sim.comm )
            sim.ptcl[1].track( sim.comm )
    else:
        # Load the fields and particles from the latest checkpoint file
        restart_from_checkpoint( sim )

    # Configure the moving window
    sim.set_moving_window( v=v_window )

    # Add a field diagnostic
    sim.diags = [ FieldDiagnostic( diag_period, sim.fld, comm=sim.comm ),
                ParticleDiagnostic( diag_period, {"electrons" : sim.ptcl[0]},
                                select={"uz" : [1., None ]}, comm=sim.comm ),
                ParticleDiagnostic( diag_period, {"nitrogen" : sim.ptcl[1]},
                                select={"uz" : [1., None ]}, comm=sim.comm )]
    # Add checkpoints
    if save_checkpoints:
        set_periodic_checkpoint( sim, checkpoint_period )

    ### Run the simulation
    sim.step( N_step )
    print('')
