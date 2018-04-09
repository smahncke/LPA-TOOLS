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
	plt.figure(figsize=size)
	for i in range(0,len(input_data)):
		try:
			plt.plot(z,input_data[i])
		except ValueError:
			print("")
			print("**NOTE: The "+str(i+1)+". dataset has different x/y dimensions - Plotting impossible")
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	title(plot_title)

	plt.show()


