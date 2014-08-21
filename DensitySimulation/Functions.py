from __future__ import division

import math as m

import numpy as np
from numpy.lib.scimath import sqrt
from scipy import constants
from scipy.integrate import dblquad

import LinAl as LA


""" This class contains the functions we have in analytical form. Contains non-static methods only because
of dependency on the initial parameters.
"""
class AnalyticalFunctions(object):
	def __init__(self, iP, dC):
		# starting points
		self.iP = iP
		self.dC = dC

	""" unperturbed surface density sigma0. uses sigmaStar. definition page 964, 23a.
	"""
	def sigma0(self,r):
		return self.dC.sigmaStar*(self.iP.rStar/float(r))**self.iP.p

	""" sound speed in the disk. definition page 964 before 23a in the bulk block of text.
	"""
	def a0(self,r):
		return self.iP.a0Star*float(r)**(-self.iP.q/2)

	""" the capital sigma. definition page 961, 1c
	"""
	def Sigma(self,r):
		return (2*m.pi*constants.G*self.sigma0(float(r)))/(self.a0(float(r))**2)


""" this class represents everything we have in a discrete manner: discretized analytic functions and discrete
functions. Besides that, a - also public available - discretizer is available based on the given initial parameters.
A discrete function will always return the whole set of values in a numpy array (due to performance reasons not
a list). thus, they need to be created and initialized what will not be done when created because it is a
time and computing intense action.
Besides that, there's a method for return the j-th component for a given discretized function - just for convenience.
"""
class DiscreteFunctions(object):
	def __init__(self, analyticFcts):
		# simple assignments for convenience
		self.iP = analyticFcts.iP
		self.dC = analyticFcts.dC
		self.analyticFunctions = analyticFcts
		# create a local version of the linear algebra class
		# by manipulating this field, the user could use a different
		# linAlg implementation because there is nothing like an defined
		# interface
		self.linAlg = LA.LinearAlgebraFunctions(self.iP, self. dC)
		self.linAlg.initLDMOne()
		self.linAlg.initLDMTwo()
		# the former analytic functions: will be numpy arrays
		self.a0Discrete = None
		self.SigmaDiscrete = None
		self.sigma0Discrete = None
		# the discrete functions right from the beginning
		self.Omega = None
		self.DOmega=None
		self.DDOmega=None
		self.kappa = None
		self.Dkappa=None
		self.initA0Discrete()
		self.initSigma0Discrete()
		self.initOmega()
		self.initKappa()
		self.initSigmaDiscrete()

	# =============== discretizing function ===============

	""" for the given function, the set of values in respect to the
	grid will be returned.
	"""
	def discretizer(self, function):
		res = np.array([])
		# iterates over the rValues and saves the
		# function value.
		for cellR in self.dC.radialCells.rValues:
			res = np.append(res, function(cellR))
		return res

	# ================ component function =================

	""" for a given set of discrete function values, the j-th
	component will be returned. for convenience purposes.
	There is no exception handling regarding index boundaries,
	so watch out.
	besides that, here the papers convention of 1-indexing is used.
	"""
	def returnComponent(self, setOfDiscreteFunctionValues, jComponent):
		return setOfDiscreteFunctionValues[jComponent-1]

	# ================ discrete functions =================
	""" calculates the set of values for the Omega function. It's the
	angular velocity profile.
	definition page 962, 5a and page 972, A1-A5.
	"""
	def calcOmega(self):
	#returns discretized capital omega as one dimensional numpy array. See appendix A.

		def omega_continuous(r):
			# create continuous function of r to feed through discretizer.
			#Actually equal to Omega^2

			def omega_star(r):
				#Equation A2 solved for omega star. Modified to use softened gravity.
				return constants.G*self.iP.mStar/(r+self.iP.eta)**2

			def omega_disk(r):
				#Equation A3 and A5. Softened gravity used throughout.

				def func(phi,x,r):
					#integrand from A3 and A5.
					return (constants.G*self.analyticFunctions.sigma0(r*x)*x)/(1+(x**2)-2*x*m.cos(phi)+self.iP.eta**2)**-.5

				def upper_bound(x):
					return float(2*m.pi)
				def lower_bound(x):
					return float(0)

				return r*(dblquad(func,self.iP.rStar/r,self.iP.rDisk/r,lower_bound,upper_bound,args=(r,))[0])

			def omega_pressure(r):
				#Equation A4
				A=(self.analyticFunctions.a0(r)**2)/self.analyticFunctions.sigma0(r) #a_0^2/sigma_0 term
				B=(-self.iP.p*self.dC.sigmaStar*self.iP.rStar**self.iP.p)/r**(self.iP.p+1) #d sigma_0/dr term
				return A*B


			return omega_star(r)/(r+self.iP.eta)+omega_disk(r)/(r+self.iP.eta)+omega_pressure(r)/r

		return sqrt(-1*self.linAlg.firstDerivative(self.discretizer(omega_continuous)))

	"""
	Calculates the first derivative of capital omega using numerical differentiation method in ARS. Returns vector with
	components equaling value at each cell.
	"""
	def calcDOmega(self):
		return self.linAlg.firstDerivative(self.Omega)

	"""
	Calculates the second derivative of capital omega using numerical differentiation method in ARS. Returns vector with
	components equaling value at each cell.
	"""
	def calcDDOmega(self):
		return self.linAlg.secondDerivative(self.Omega)

	""" calculates the set of values for the kappa function. uses the
	derivative matrices to form d/dr.
	definition page 962, 5b.a
	"""
	def calcKappa(self):
		prefactor = self.dC.radialCells.rValues**(-4.0)
		toDerive = np.asmatrix(((self.dC.radialCells.rValues**2)*self.Omega)**2).T
		derived = (np.asarray((self.linAlg.logDerivationMatrixOne * toDerive).T)).reshape([len(prefactor)])
		squared = (prefactor * derived)
		return sqrt(squared)
		# TODO: returns interesting set of values, only boundary conditions differ
		# even tough 1/r^4!
		# TODO: can kappa be element of complex number?
		# TODO: bring into line with current derivation method. Use DOmega.

	"""
	Calculates the first derivative of kappa analytically. As kappa is dependent on omega, which is discrete, the result is also discrete. Returns vector with
	components equaling value at each cell.
	"""
	def calcDKappa(self):
		def calcDKappaAnalytic(r):
			return .5*(4*self.Omega(r)**2+2*r*self.Omega(r)*self.DOmega(r))**-.5*(10*self.Omega(r)*self.DOmega(r)+2*r*self.DOmega(r)**2+2*r*self.Omega(r)*self.DDOmega(r))

		return self.discretizer(calcDKappaAnalytic)

	# ============== init discretized fields ==============
	""" call this function in order to initialize the discrete values for
	a0. if already initialized, an error is raised.
	uses the discretizer acting upon the analytic function.
	"""
	def initA0Discrete(self):
		if self.a0Discrete == None:
			self.a0Discrete = self.discretizer(self.analyticFunctions.a0)
		else:
			raise BaseException("a0 already initialized")

	""" call this function in order to initialize the discrete values for
	capital Sigma. if already initialized, an error is raised.
	uses the discretizer acting upon the analytic function.
	"""
	def initSigmaDiscrete(self):
		if self.SigmaDiscrete == None:
			self.SigmaDiscrete = self.discretizer(self.analyticFunctions.Sigma)
		else:
			raise BaseException("Sigma already initialized")

	""" call this function in order to initialize the discrete values for
	 sigma0. if already initialized, an error is raised.
	uses the discretizer acting upon the analytic function.
	"""
	def initSigma0Discrete(self):
		if self.sigma0Discrete == None:
			self.sigma0Discrete = self.discretizer(self.analyticFunctions.sigma0)
		else:
			raise BaseException("sigma0 already initialized")

	""" call this function in order to initialize the discrete values for
	capital Omega. if already initialized, an error is raised.
	uses the numerical integral method given in the paper.
	"""
	def initOmega(self):
		if self.Omega == None:
			self.Omega = self.calcOmega()
		else:
			raise BaseException("Omega already initialized")

	"""call this function in order to initialize the discrete values for
	the first derivative of capital Omega. if already initialized, an error is raised.
	uses the numerical integral method given in the paper."""

	def initDOmega(self):
		if self.DOmega == None:
			self.DOmega = self.calcDOmega()
		else:
			raise BaseException("DOmega already initialized")

	"""call this function in order to initialize the discrete values for
	the second derivative of capital Omega. if already initialized, an error is raised.
	uses the numerical integral method given in the paper."""

	def initDDOmega(self):
		if self.DDOmega == None:
			self.DDOmega = self.calcDDOmega()
		else:
			raise BaseException("DDOmega already initialized")

	""" call this function in order to initialize the discrete values for
	kappa. if already initialized, an error is raised.
	matrix operations are used to perform the derivation.
	"""
	def initKappa(self):
		if self.kappa == None:
			self.kappa = self.calcKappa()
		else:
			raise BaseException("kappa already initialized")

	""" call this function in order to initialize the discrete values for the first derivative of
	kappa. if already initialized, an error is raised.
	derivation is analytical, but dependence on omega means values are discrete.
	"""
	def initDKappa(self):
		if self.DKappa == None:
			self.DKappa = self.calcDKappa()
		else:
			raise BaseException("DKappa already initialized")





# ======================== TESTING AREA ========================

"""import Params as P
iP = P.InitialParameters(number = 10)
dC = P.DerivedConstants(iP)
anF = AnalyticalFunctions(iP, dC)
diF = DiscreteFunctions(anF)

diF.initOmega()
diF.initKappa()

print diF.kappa
print "hi"

"""