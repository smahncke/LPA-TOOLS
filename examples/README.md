# PDP-EXAMPLE

# Short description:

This example plots the plasma density profile of a hydrogen plasma with a sin^2-like, 500e-6 meter long up- and downramp and a 
constant, 4mm long plateau. The nitrogen density has the shape of a gauss curve with FWHM (full with at half max) of 300e-6 m. 

# How to use it:

NOTE: If you want to change parameters like the resolution of the simulation outputs, you have to do it in the 'FBPIC_example.py'
      file (see https://github.com/fbpic/fbpic for further information). Otherwise you can easily use the already existing 'diags' folder.
      Make sure Python is correctly installed on your computer.
      
1.  Download the 'FBPIC_example.py' and the 'pdp_example.py' file. If you don't need to have a higher simulation resolution of the
    Particle-In-Cell code FBPIC you can also just download the 'diags' folder and jump to the third step.
2.  Run 'FBPIC_example.py' the terminal by typing
`python FBPIC_example.py` 
. You have to install FBPIC first, see the documentation.
3.  If you haven't done yet, install openPMD_viewer (see https://github.com/openPMD/openPMD-viewer)
4.  Now run 'pdp_example.py' in the terminal by typing `phyton pdp_example.py` and follow the instructions.

