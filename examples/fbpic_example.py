"""
This is a typical input script that runs a simulation of
laser-wakefield acceleration using FBPIC.

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
where fbpic_object is any of the objects or function of FBPIC.
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
Nz = 400         # Number of gridpoints along z 
zmax = 30.e-6    # Right end of the simulation box (meters)
zmin = -60.e-6   # Left end of the simulation box (meters)
Nr = 100          # Number of gridpoints along r
rmax = 100.e-6    # Length of the box along r (meters)
Nm = 2           # Number of modes used

# The simulation timestep
dt = (zmax-zmin)/Nz/c   # Timestep (seconds)
N_step = int(5.1e-3/(dt*c))    # Number of iterations to perform #ORIGINAL int(1.2e-3/(dt*c))


# Order of accuracy of the spectral, Maxwell (PSATD) solver.
# -1 correspond to infinite order, i.e. wave propagation is perfectly
# dispersion-free in all directions. This is adviced for single GPU/CPU
# simulations. For multi GPU/CPU simulations, choose n_order > 4
# (and multiple of 2). A large n_order leads to more overhead in MPI
# communications, but also to a more accurate dispersion for waves.
# (Typically, n_order = 32 gives accurate physical results)
n_order = -1

# The gas
dope_element = 'N' # Element the H2 will be doped with (i.e. 'N' or 'He')
n_init = 0.5e18*1.e6 # Gas density of the doped gas
am_H2 = 99. # Amount of H2 particles in percent
am_dope_gas = 1.  # Amount of dope gas particles in percent (sum of am_H2 and am_dope_gas must be 100 %)
pre_ion_level = 5. # Pre-ionization level of the dope gas atoms
q_dope_gas = e # Charge of the ionized particles (standard value is the electron charge e)

n_H2 = (am_H2/100)*n_init #Gas density of H2
n_dope_gas = (am_dope_gas/100)*n_init #Gas density of dope gas

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

# The density profile
ramp_start = 30.e-6
ramp_length = 500.e-6
plateau = 4000.e-6
def dens_func( z, r ) :
    """Returns relative density at position z and r"""
    # Allocate relative density
    n = np.ones_like(z)
    n = np.where( z<ramp_start+ramp_length, np.sin(
                 (0.5*np.pi*(z-ramp_start))/ramp_length)**2, n )
    n = np.where( z>ramp_start+ramp_length+plateau, np.sin(
                 (0.5*np.pi*(z-ramp_start))/ramp_length)**2, n )
    # Supress density before and after the ramp
    n = np.where( z<ramp_start, 0., n )
    n = np.where( z>ramp_start+2*ramp_length+plateau, 0., n)
    return(n)

# Specifies of the dope gas

if dope_element == 'N': # Gives the mass of the according dope gas
    m_dope_element = 14.*m_p
elif dope_element == 'He':
    m_dope_element = 4.*m_p
elif dope_element == 'Ar':
    m_dope_element = 40.*m_p
elif dope_element == 'Kr':
    m_dope_element = 84.*m_p
else:
    print('ERROR: Element '+ dope_element + ' is not defined at the moment. Please choose argon, krypton, nitrogen or helium ')

if dope_element == 'N': # For N2 (or generally for no noble gases) the amount of ionized electrons will be multiplied by 2
    molecular_factor = 2.
else:
    molecular_factor = 1.

# The plasma

n_e_H2 = n_H2*2. #Plasma density of the H2 (without the pre-ionized electrons)
n_dope_gas_pre_ion = n_dope_gas*pre_ion_level*molecular_factor #Electron density of the pre-ionized dope gas (depends on the pre-ionization level)

n_e = n_e_H2+n_dope_gas_pre_ion #Plasma density of H2 including the pre-ionized electrons


# ---------------------------
# Carrying out the simulation
# ---------------------------

# NB: The code below is only executed when running the script,
# (`python -i lpa_sim.py`), but not when importing it (`import lpa_sim`).
if __name__ == '__main__':


    # Initialize the simulation object
    sim = Simulation( Nz, zmax, Nr, rmax, Nm, dt,
        p_zmin, p_zmax, p_rmin, p_rmax, p_nz, p_nr, p_nt, n_e,
        dens_func=dens_func, zmin=zmin, boundaries='open',
        n_order=n_order, use_cuda=use_cuda )
    
    # Add the N atoms
    p_zmin, p_zmax, Npz = adapt_to_grid( sim.fld.interp[0].z,
                                        p_zmin, p_zmax, N_nz )
    p_rmin, p_rmax, Npr = adapt_to_grid( sim.fld.interp[0].r,
                                        p_rmin, p_rmax, N_nr )
    sim.ptcl.append(
        Particles(q_dope_gas, m_dope_element, n_dope_gas, Npz, p_zmin,
            p_zmax, Npr, p_rmin, p_rmax,
            N_nt, dt, use_cuda=use_cuda, dens_func=dens_func,
            grid_shape=sim.fld.interp[0].Ez.shape,
            continuous_injection=True ) )
            
    sim.ptcl[1].make_ionizable(element=dope_element, level_start=pre_ion_level,
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
                ParticleDiagnostic( diag_period, {"dope element (" + dope_element + ")" : sim.ptcl[1]},
                                select={"uz" : [1., None ]}, comm=sim.comm )]
    # Add checkpoints
    if save_checkpoints:
        set_periodic_checkpoint( sim, checkpoint_period )

    ### Run the simulation
    sim.step( N_step )
    print('')
