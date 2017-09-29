# Plasma density plotter (pdp)

# DESCRIPTION:

PDP is an interactive script, which plots plasma density profile for ionization-induced trapping laser plasma accelerators (see an example of a plot below).

NOTE: This code is just for experimental use at the moment. It's not finished yet.

# REQUIREMENTS:

- A python distribution should be installed to run the code (i.e. Anaconda, see https://docs.anaconda.com/)
- The Particle-In-Cell code FBPIC should be installed to create the (laser) diagnostic files. For further information about the 
  installation of FBPIC click on the following link: https://github.com/fbpic/fbpic
- You'll also need the tool 'openPMD-viewer' for the laser amplitude analytics. For further information click on the following
  link: https://github.com/openPMD/openPMD-viewer


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
3. The next step depends on the version of the pdp-script:
  a)  tools/pdp_main.py: This is the original version of the plotting script, you have to set some of the global parameters (like    
      initial amount of nitrogen gas or pre-ionization level) manually. To do so, open the python file in a text editor and set     
      the values of these parameters in a way you want them to be. 
  b)  pdp_GUI.py: This is the newer, more interactive version of this script. It does the same like the original one, but the 
      paramaters can be changed directly in the terminal while the script is running. 
4. Start the script in a terminal and follow the instructions.
5. That's it.
      
# TOOLS:

In the 'tools' directory there are some more pdp-scripts. Read the README.md-file in it for more information.

# KNOWN BUGS:

- At the moment the whole plot is shifted a bit to lower values of z [mm], that's why the hydrogen curve has the short upramp
  for z > 5mm
- The conversion of the laser amplitude into a energy isn't correct, but this bug will be fixed in the near future

# EXAMPLE

![alt Example](https://github.com/smahncke/pdp/blob/master/first_draft.png?raw=true)

# OLDER VERSIONS

The older versions of this script will be saved in the 'source' folder. 

Feel free to contact me if you need further information or if you got some issues: sebastian.mahncke@desy.de


