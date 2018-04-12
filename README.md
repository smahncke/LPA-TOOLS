# Laser-plasma-acceleration tools (LPA-TOOLS)

## Short description: 

This repository includes a bunch of simple tools for simulations of laser-plasma-acceleration, especially ionization-induced trapping in laser-plasma-accelerators.

## Features:

This is a short overview of the features that are included in this repository. More detailed README files of different features will be added in the correlating directory soon.

### Ionization induced trapping

For ionization induced trapping it's important to know the degree of ionization, after the laser pulse propagated through the gas.
The degree of ionization depends on different parameters like the element that is used for doping the plasma target and a bunch of
laser probabilities. With the LPA-Tool, you can calculate the 
- Ionization probability  of the dope gas within the laser pulse
- The degree of ionization of the dope gas at the end and at an arbitary position of the laser pulse

The results can be used to estimate a minimal needed laser intensity/ a0 to run a laser plasma accelerator. You can also:
- calculate and plot the wakefield 
- calculate the so-called 'trapping condition', which is essential for ionization-induced trapping. It can be used to calculate analytically the trapped charge/bunch charge (this feature will be added soon, too).

#### Example

![alt example_plot](example_plot.png)

## How to install and use the script:

#### Requirements

All tools are python scripts, so you should've installed a python distribution, for example Anaconda, see https://www.anaconda.com/download/

#### How to install the code 

After installing Anaconda, you can install LPA-TOOLS by making the following steps:

1. Open a terminal and browse to the directory you want to install this script.
2. Clone this repository by typing: 
	```
	git clone https://github.com/smahncke/LPA-TOOLS
	```	

#### How to use the code

1. You can use the LPA-TOOLS_example.py script for example or modify it with a text editor. Take a look at the documentation for the syntax.
2. Open a terminal and browse to the directory you installed LPA-TOOLS
3. Run the script by typing (replace 'script_name' by the actual name of the script):
	```
	python script_name.py
	```




## Contact:

Feel free to contact me if you need further information or if you got some issues: sebastian.mahncke@desy.de

Thank you for using LPA-TOOLS. 

