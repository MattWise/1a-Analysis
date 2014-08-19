from __future__ import division

import math as m
import Params as P
from scipy import constants
from scipy.integrate import dblquad
import numpy as np

""" This class contains the functions we have in analytical form. Contains non-static methods only because
of dependency on the initial parameters.
"""
class AnalyticalFunctions(object):
	def __init__(self, iP, dC):
		self.iP = iP
		self.dC = dC

	# unperturbed disk surface density
	def sigma0(self,r):
		return self.dC.sigmaStar*(self.iP.rStar/float(r))**self.iP.p

	# sound speed
	def a0(self,r):
		return self.iP.a0Star*float(r)**(-self.iP.q/2)

	# big sigma
	def sigma(self,r):
		return (2*m.pi*self.iP.G*self.sigma0(float(r)))/(self.a0(float(r))**2)


""" This class contains the functions we have in discrete appearance. Contains non-static methods only because
of dependency on the initial parameters and uses discretized analytical methods by using the (static) functions
in the related classes.
"""
class DiscreteFunctions(object):
	def __init__(self, analyticFcts):
		self.iP = analyticFcts.iP
		self.dC = analyticFcts.dC
		self.analyticFunctions = analyticFcts

	# kappa
	def kappa(rJ):
		pass

	def unperturbed_rotational_velocity(self):
	#returns discretized capital omega as one dimensional numpy array. See appendix A.

		def omega_continuous(r):
			# create continuous function of r to feed through discretizer.
			#Actually equal to r*Omega^2

			def omega_star(r):
				#Equation A2 solved for omega star. Modified to use softened gravity.
				return (constants.G*self.iP.MStar/(r+self.iP.softening_parameter)**3)**.5

			def omega_disk(r):
				#Equation A3 and A5. Softened gravity used throughout.

				def func(phi,x,r):
					#integrand from A3 and A5.
					return (constants.G*self.analyticFunctions.sigma_0(r*x)*x)/(1+(x**2)-2*x*m.cos(phi)+self.iP.softening_parameter**2)**-.5

				def upper_bound(x):
					return float(2*m.pi)
				def lower_bound(x):
					return float(0)

				return r*(dblquad(func,self.iP.RStar/r,self.iP.RDisk/r,lower_bound,upper_bound,args=(r,))[0])

			def omega_pressure(r):
				#Equation A4
				A=(self.analyticFunctions.a0(r)**2)/self.analyticFunctions.sigma_0(r) #a_0^2/sigma_0 term
				B=(-self.iP.p*self.dC.sigmaStar*self.iP.RStar**self.iP.p)/r**(self.iP.p+1) #d sigma_0/dr term
				C=(A*B/r)**.5 #solve for omega_P
				return C


			return omega_star(r)+omega_disk(r)+omega_pressure(r)

		def divide_by_r(array):
			#divides through by r, yielding discretized omega^2
			for index,value in enumerate(self.dC.radial_cells.RList):
				array[index]=array[index]/value

		return np.sqrt(divide_by_r(-1*first_derivative_place_holder(self.discretizer(omega_continuous))))

