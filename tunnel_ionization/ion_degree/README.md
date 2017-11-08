# Laser-plasma-acceleration tools (LPA-TOOLS)

## Tunnel ionization degree (TID)

### DESCRIPTION:

TID is a python script that plots the degree of ionization of a gas assuming tunnel ionization probability in dependence of typical parameters used for laser plasma acceleration, i.e. laser amplitude, laser profile etc.

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

- A plot that shows the electric field of the gaussian laser pulse (dotted blue line), the laser envelope (orange line), the tunnel ionization probability (red line) which depends on the laser pulse and the potential of the ionization level, and the degree of ionization (green line). The laser propagates from right (positive z values) to the left. 

- The max ionization degree in percent
- The max ionization probability in percent
- Laser amplitude threshold a(z) at which the ionization starts (ionization probability >= 1%)

### KNOWN BUGS/ISSUES:

(There are no known issues at the moment. Feel free to contact me (sebastian.mahncke@desy.de) if you found some bugs)

### EXAMPLES:
(Click on the image for a bigger view or browse into the examples folder)
1. Example of nitrogen (level N4+ -> N5+):
![alt Example](https://github.com/smahncke/LPA-TOOLS/blob/master/tunnel_ionization/ion_degree/examples/Nitrogen/N_5.png?raw=true)

2. Example of argon (level Ar9+ -> Ar10+):
![alt Example](https://github.com/smahncke/LPA-TOOLS/blob/master/tunnel_ionization/ion_degree/examples/Argon/Ar_10.png?raw=true)

### OTHER TOOLS:

1. Automated: An automated version that creates .png-files for ALL possible ionization levels automatically.


Feel free to contact me if you need further information or if you got some issues: sebastian.mahncke@desy.de


