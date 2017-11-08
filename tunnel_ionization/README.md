# Laser-plasma-acceleration tools (LPA-TOOLS)

## Tunnel ionization probability (TIP)

### DESCRIPTION:

TIP is a python script that plots the tunnel ionization probability due to laser ionization for different gases in 
dependence of typical parameters used for laser plasma acceleration, i.e. laser amplitude.

### REQUIREMENTS:

- A python distribution should be installed to run the code (i.e. Anaconda, see https://docs.anaconda.com/)

### FEATURES:

- Nitrogen, Argon and Krypton are choosable

### HOW TO USE:

1. Make sure the requirements are fulfilled.
2. Open the python file in a text editor and set the initial values in the header of the script. The corresponding
   parameters depend on the script you want to start.
3. Start the script in a terminal by typing 'python script_name.py' (replace script_name.py by the actual name of the
   script you want to start).
4. That's it.
      
### OUTPUTS:

Depends on the script:

- 'energy_prob_2D.py': 2D plot that gives the tunnel ionization probability in dependence of the initial energy of the ionized electrons after ionization. epsilon = 0 f.e. means that the electron will be born at rest.
- 'laser_prob_2D.py': 2D plot that gives the tunnel ionization probability in dependence of the laser amplitude.
- 'energy_laser_prob_3D': Combines both of the scripts above -> 3D plot that gives the ionization probability in dependence of the initial electron energy AND the laser amplitude. 

### KNOWN BUGS/ISSUES:

(There are no known issues at the moment. Feel free to contact me (sebastian.mahncke@desy.de) if you found some bugs)

### EXAMPLES:
(Click on the image)
1. Example of 'laser_prob_2D':
![alt Example](https://github.com/smahncke/LPA-TOOLS/blob/master/tunnel_ionization/ion_probability/examples/example_2D_Argon.png?raw=true)

2. Example of 'energy_laser_prob_3D':
![alt Example](https://github.com/smahncke/LPA-TOOLS/blob/master/tunnel_ionization/ion_probability/examples/example_3D_Argon.png?raw=true)

### OTHER TOOLS:

1. Automated: An automated version of the 3D plot. Creates .png-files for ALL possible ionization levels automatically.


Feel free to contact me if you need further information or if you got some issues: sebastian.mahncke@desy.de


