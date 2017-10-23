# Plasma density analysis (PDA)

# DESCRIPTION:

PDA is a simple python script that plots gas density profiles for ionization-induced trapping laser plasma accelerators (see an example of a plot below).

# REQUIREMENTS:

- A python distribution should be installed to run the code (i.e. Anaconda, see https://docs.anaconda.com/)

# FEATURES:

- The initial amount of N2 can be set arbitrarily. 
- The density profile of the nitrogen corresponds to a gauss profile. The FWHM (full-width at half max) of this 
  gauss profile can be set arbitrarily
- The nitrogen gas can be chosen pre-ionized or non pre-ionized
- Some parameters like the upramp length and the length of the hydrogen density profile can be set manually


# HOW TO USE:

1. Make sure the requirements are fulfilled.
2. Create a bunch of .hdf5 diagnostic files by running the Particle-In-Cell code FBPIC (see above). Make sure that the 'diags' 
   folder is in the same directory like the pdp script.
3. Open the python file in a text editor and set the initial values of the FWHM, final electron density and nitrogen amount in a way you want them to be. 
4. Start the script in a terminal by typing 'python pda_main.py'.
5. That's it.
      
# OUTPUTS:

There are some output files that will be saved in the 'output' folder:

1. The plot of the gas density profiles
2. The arrays of the gas densities (nitrogen density, hydrogen density and the total density)

# KNOWN BUGS/ISSUES:

(There are no known issues at the moment. Feel free to contact me (sebastian.mahncke@desy.de) if you found some bugs)

# EXAMPLE

![alt Example](https://github.com/smahncke/pdp/blob/master/first_draft.png?raw=true)

# OLDER VERSIONS

The older versions of this script will be saved in the 'source' folder. 

Feel free to contact me if you need further information or if you got some issues: sebastian.mahncke@desy.de


