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
		return np.complex128(self.dC.sigmaStar*(self.iP.rStar/float(r))**self.iP.p)

	"""first derivative of surface density sigma0. Calculated analytically.
	"""
	def Dsigma0(self,r):
		return np.complex128((-self.iP.p*self.dC.sigmaStar*self.iP.rStar**self.iP.p)/(r**(self.iP.p+1)))

	""" sound speed in the disk. definition page 964 before 23a in the bulk block of text.
	"""
	def a0(self,r):
		return np.complex128(self.iP.a0Star*float(r)**(-self.iP.q/2))

	""" the capital sigma. definition page 961, 1c
	"""
	def Sigma(self,r):
		return np.complex128((2*m.pi*constants.G*self.sigma0(float(r)))/(self.a0(float(r))**2))


""" this class represents everything we have in a discrete manner: discretized analytic functions and discrete
functions. Besides that, a - also public available - discretizer is available based on the given initial parameters.
A discrete function will always return the whole set of values in a numpy array (due to performance reasons not
a list). thus, they need to be created and initialized what will not be done when created because it is a
time and computing intense action. there is a init() method available which calls all the init methods in order to have
a user-convenient way. but nevertheless, all the init methods can be called separately if needed.
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
		# the former analytic functions: will be numpy arrays
		self.a0Discrete = None
		self.SigmaDiscrete = None
		self.sigma0Discrete = None
		self.Dsigma0Discrete=None
		# the discrete functions right from the beginning
		self.Omega = None
		self.DOmega=None
		self.DDOmega=None
		self.kappa = None
		self.Dkappa = None

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

	# ================ discrete functions =================
	""" calculates the set of values for the Omega function. It's the
	angular velocity profile.
	definition page 962, 5a and page 972, A1-A5.
	"""
	def calcOmega(self):
	#returns discretized capital omega as one dimensional numpy array. See appendix A.

		def psi_grav(r):
			# create continuous function of r to feed through discretizer.
			# potential due to star and disk divided by r.

			def psi_star(r):
				#Equation A2 solved for omega star. Modified to use softened gravity.
				return -constants.G*self.iP.mStar/(r+self.iP.eta)

			def psi_disk(r):
				#Equation A3 and A5. Softened gravity used throughout.

				def func(phi,x,r):
					#integrand from A3 and A5.
					return (constants.G*self.analyticFunctions.sigma0(r*x)*x)/(1+(x**2)-2*x*m.cos(phi)+self.iP.eta**2)**-.5

				def upper_bound(x):
					return float(2*m.pi)
				def lower_bound(x):
					return float(0)

				return r*-(dblquad(func,self.iP.rStar/r,self.iP.rDisk/r,lower_bound,upper_bound,args=(r,))[0])

			return np.complex128((psi_star(r)+psi_disk(r)/r))

		def psi_pressure(r):
			#Equation A4
			#Actually first derivative with respect to r of psi pressure
			A=(self.analyticFunctions.a0(r)**2)/self.analyticFunctions.sigma0(r) #a_0^2/sigma_0 term
			B=(-self.iP.p*self.dC.sigmaStar*self.iP.rStar**self.iP.p)/r**(self.iP.p+1) #d sigma_0/dr term
			return A*B

		psi_grav_discrete=self.linAlg.firstDerivative(self.discretizer(psi_grav))
		psi_pressure_discrete=self.discretizer(psi_pressure)

		return sqrt(psi_grav_discrete+psi_pressure_discrete)

	""" Calculates the first derivative of capital omega using numerical differentiation method in ARS. Returns vector with
	components equaling value at each cell.
	"""
	def calcDOmega(self):
		return self.linAlg.firstDerivative(self.Omega)

	""" Calculates the second derivative of capital omega using numerical differentiation method in ARS. Returns vector with
	components equaling value at each cell.
	"""
	def calcDDOmega(self):
		return self.linAlg.secondDerivative(self.Omega)

	""" calculates the set of values for the kappa function. uses the
	derivative matrices to form d/dr.
	definition page 962, 5b.a
	"""
	def calcKappa(self):
		array=np.zeros(self.iP.number, dtype=np.complex128)
		for index,r in enumerate(self.dC.radialCells.rValues):
			array[index]=np.complex128(4*self.Omega[index]**2+2*r*self.Omega[index]*self.DOmega[index])
		return sqrt(array)

	"""	Calculates the first derivative of kappa analytically. As kappa is dependent on omega, which is discrete, the result is also discrete. Returns vector with
	components equaling value at each cell.
	"""
	def calcDKappa(self):
		array1=(.5*self.kappa)**-1
		array2=np.zeros(self.iP.number,dtype=np.complex128)
		for index,r in enumerate(self.dC.radialCells.rValues):
			array2[index]=np.complex128(10*self.Omega[index]*self.DOmega[index]+2*r*self.DOmega[index]**2+2*r*self.Omega[index]*self.DDOmega[index])
		return array1*array2

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

	"""call this function in order to initialize the discrete values for
	 the first derivative of sigma0. if already initialized, an error is raised.
	uses the discretizer acting upon the analytic function.
	"""
	def initDSigma0Discrete(self):
		if self.Dsigma0Discrete == None:
			self.Dsigma0Discrete = self.discretizer(self.analyticFunctions.Dsigma0)
		else:
			raise BaseException("Dsigma0 already initialized")

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
	def initDkappa(self):
		if self.Dkappa == None:
			self.Dkappa = self.calcDKappa()
		else:
			raise BaseException("Dkappa already initialized")

	"""call this function to initialize first derivative matrix. Must be initialized before DOmega and Dkappa."""

	def initLDMOne(self):
		if self.linAlg.logDerivationMatrixOne==None:
			self.linAlg.initLDMOne()
		else:
			raise BaseException("logDerivationMatrixOne already initialized")

	"""call this function to initialize second derivative matrix. Must be initialized before DDOmega."""

	def initLMDTwo(self):
		if self.linAlg.logDerivationMatrixTwo==None:
			self.linAlg.initLDMTwo()
		else:
			raise BaseException("logDerivationMatrixTwo already initialized")

	"""call this function to initialize the j integration matrix. Must be initialized before proceeding to linear algebra."""

	def initJMatrix(self):
		if self.linAlg.JMatrix==None:
			self.linAlg.initJMatrix()
		else:
			raise BaseException("JMatrix already initialized")


	# ======================== Init Function ========================
	""" call this function in order to initialize the discrete values for
	all functions (i.e.: kappa, Omega, sigma0, sigma, a0).
	"""
	def init(self):
		self.initLDMOne()
		self.initLMDTwo()
		self.initSigma0Discrete()
		self.initDSigma0Discrete()
		self.initA0Discrete()
		self.initOmega()
		self.initDOmega()
		self.initDDOmega()
		self.initKappa()
		self.initDkappa()
		self.initSigmaDiscrete()
		self.initJMatrix()


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