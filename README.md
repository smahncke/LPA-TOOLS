# pdp
A plasma density profile plotter for ionization-induced trapping in laser plasma accelerator
--------------------------------------------------------------------------------------------
NOTE: Requires a bunch of laser diagnostic files created by a 
spectral, quasi-3D Particle-In-Cell code (f.e. FBPic, see https://github.com/fbpic/fbpic)
--------------------------------------------------------------------------------------------

# FEATURES:

- The initial amount of N2 can be set arbitrarily. 
- The density profile of the nitrogen corresponds to a gauss profile. The FWHM (full-width at half max) of this 
  gauss profile can be set arbitrarily
- The nitrogen gas can be chosen pre-ionized or non pre-ionized
- Some parameters like the upramp length and the length of the hydrogen density profile can be set manually


# HOW TO USE:

1. Create a bunch of .hdf5 diagnostic files with f.e. FBPic (see the link above)
2. Open the python file in a text editor and set the values of the parameters for the total initial gas density, the initial
  amount of nitrogen and the pre-ionization level in a way you want them to be. 
3. Start the script in a terminal, the 'diags/hdf5' folder created by fbpic should be in the same folder like the pdp-script
4. That's it.

NOTE: In the near future a more interactive script will be released in which some parameters can be changed directly in the 
      terminal, so the second step above is no longer necessary. 
      
# TOOLS:

In the 'tools' directory there are some more pdp-scripts. Read the README.md-file in it for more information.

# KNOWN BUGS:

- At the moment the whole plot is shifted a bit to lower values of z [mm], that's why the hydrogen curve has the short upramp
  for z > 5mm
- The conversion of the laser amplitude into a energy isn't correct, but this bug will be fixed in the near future

# Example

[[ https://github.com/smahncke/pdp/blob/master/first_draft.png?raw=true|alt=Example]]
