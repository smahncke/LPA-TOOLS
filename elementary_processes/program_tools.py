import numpy as np
from scipy.integrate import odeint
from pylab import *
import scipy.constants as const
from numba import double
from numba.decorators import jit, autojit
import scipy.constants as const
import scipy.integrate as integrate
import sys
sys.path.append('../')


def plotter(z,input_data,size = (15,10),x_label = "", y_label = "", plot_title = ""):
	"""
	Plots the given input data

	-INPUT:	-z [array]:				z axis of the plots
			-input_data [array]:	Each entry of this array will be plotted and has to have the same dimensions like 'z'
			-size [2d array,
					optional]:		Size of the plot [width, height]
			-x_label [string,
					optional]:		Label of the x axis
			-y_label [string,
					optional]:		Label of the y axis
			-plot_title [string,
					optional]:		Plot of the label

	-RETURN:- Plots [Images]
	"""

	#Initialize the plot with the given size
	plt.figure(figsize=size)

	for i in range(0,len(input_data)):
		
		#Checks if the x and y dimensions match
		try:
			plt.plot(z,input_data[i])
		except ValueError:
			print("")
			print("**NOTE: The "+str(i+1)+". dataset has different x/y dimensions - Plotting impossible")
	
	#Labels	
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	
	#Title
	title(plot_title)

	#Open the image
	plt.show()


