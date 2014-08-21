from __future__ import division

import numpy as np
import time as t
import math as m

""" This class holds all the matrices needed for numerical derivatives and integration. Since matrix creation is
quite intense in time and power, we use a non-static class here in order to store the calculated matrices.
We use Numpy.Matrix for performance. The class needs to be initialized in order to inform the user that some
code might take some time - no exception handling here only for multiple-inits.
"""
class LinearAlgebraFunctions(object):
	def __init__(self, iP, dC):
		# basic parameters
		self.initialParameters = iP
		self.derivedConstants = dC
		# first, not initialized to save time. besides that,
		# these are the checking values if already initialized.
		# but there are no control mechanisms.
		self.logDerivationMatrixOne = None
		self.logDerivationMatrixTwo = None
		# J_ij matrix which performs an integral
		self.JMatrix = None

	# ======== init methods ========

	""" Initialization of logarithmic matrix of order 1.
	"""
	def initLDMOne(self):
		if self.logDerivationMatrixOne == None:
			self.logDerivationMatrixOne\
				= self.__logDerivationMatrix(self.initialParameters.number, 1) # use method with correct input.
		else:
			raise BaseException("logDerivationMatrix of order=1 already initialized")

	""" Initialization of logarithmic matrix of order 2.
	"""
	def initLDMTwo(self):
		if self.logDerivationMatrixTwo == None:
			self.logDerivationMatrixTwo\
				= self.__logDerivationMatrix(self.initialParameters.number, 2) # use method with correct input.
		else:
			raise BaseException("logDerivationMatrix of order=2 already initialized")

	""" Initialization of J_ij matrix.
	"""
	def initJMatrix(self):
		if self.JMatrix == None:
			self.JMatrix = self.__calcJMatrix()
		else:
			raise BaseException("JMatrix already initialized")

	""" For a given dimension size N, returns the logarithmic derivation matrices
	of order one and two: the discretized version of d/dlog(r).
	used to get d/dr and (d/dr)^2, derivation see page 973, B3a and B3b.
	"""
	def __logDerivationMatrix(self, dimension, order, innerBoundary=False, outerBoundary=False):

		#creates a differentiation matrix for numerical differentiation. "Order=1" means first derivative, "order=2" means second.


		def diagonals(array,a,b,value):
			#goes down diagonal form specified coordinates, setting each to value.
			if a>b:
				times = dimension-a-1
			else:
				times=dimension-b
			for row in range(times):
				array[a,b]=value
				a+=1
				b+=1
			return array


		def setBoundary(array,dimension,innerBoundary,outerBoundary):
			#sets the matrix boundary conditions. Defaults are from ARS paper.

			outerBoundary.reverse()

			for index,value in enumerate(innerBoundary):
				array[0,index]=value
			for index,value in enumerate(outerBoundary):
				array[dimension-1,dimension-1-index]=value
			return array

		array = np.zeros([dimension,dimension])
		prefactor = 0

		if order == 1:
			array = diagonals(array,1,0,-1)
			array = diagonals(array,1,2,1)

			#defualt boundary conditions
			if not innerBoundary:
				inner_boundary=[-3,4,-1]
			if not outerBoundary:
				outer_boundary=[1,-4,3]

			prefactor = 1/(2*np.log10(self.derivedConstants.radialCells.f))

		if order == 2:
			array = diagonals(array,1,0,1)
			array = diagonals(array,1,1,-2)
			array = diagonals(array,1,2,1)

			#defualt boundary conditions
			if not innerBoundary:
				inner_boundary = [2,-5,4,-1]
			if not outerBoundary:
				outer_boundary=[-14,-5,2]

			prefactor = 1/(np.log10(self.derivedConstants.radialCells.f))**2

		array=setBoundary(array,dimension,inner_boundary,outer_boundary)

		return prefactor * array

	def firstDerivative(self,vector):
		res=np.dot(self.logDerivationMatrixOne, vector)
		for index,value in enumerate(self.derivedConstants.radialCells.rValues):
			res[index]=res[index]/value
		return res

	def secondDerivative(self,vector):
		res=np.dot(self.logDerivationMatrixTwo, vector)-np.dot(self.logDerivationMatrixOne, vector)
		for index,value in enumerate(self.derivedConstants.radialCells.rValues):
			res[index]=res[index]/(value**2)
		return res

	""" Performs the actual J_ij matrix calculation. Since an 2d array: O(dimension^2). Actually,
	this matrix is only dependent on the r cell values.
	"""
	def __calcJMatrix(self):
		# get dimension
		dimension = self.initialParameters.number
		p = self.initialParameters.p
		printed = -5
		# get the radCells
		radCellsValues = self.derivedConstants.radialCells.rValues
		# creating the result of correct dimension filled with zeros
		result = np.zeros((dimension, dimension))
		# filling the matrix with correct values according to my calculations
		# todo: maybe parallel! would be a nice task :-)
		for i in np.arange(dimension):
			for j in np.arange(1, dimension-1):
				result[i][j] = (radCellsValues[j+1]-radCellsValues[j])/(2*radCellsValues[i])\
							   *(radCellsValues[j]/radCellsValues[i])**(2-p)\
							   + (radCellsValues[j]-radCellsValues[j-1])/(2*radCellsValues[i])\
								 *(radCellsValues[j]/radCellsValues[i])**(2-p)
			result[i][0] = (radCellsValues[1]-radCellsValues[0])/(2*radCellsValues[i])\
						   *(radCellsValues[0]/radCellsValues[i])**(2-p)
			result[i][dimension-1] = (radCellsValues[dimension-1]-radCellsValues[dimension-2])/(2*radCellsValues[i])\
						   *(radCellsValues[dimension-1]/radCellsValues[i])**(2-p)
			# some simple feedback
			percentage = round(100*(i/dimension), 0)
			if percentage%5==0 and not (percentage==printed):
				print "completed:", percentage, "%"
				printed = percentage

		return result


""" This class represents the huge matrix generated in order to solve the omega eigenvalue and
the S eigenvector.
Due to heavy computing, needs to be initialized as well.
"""
class WMatrix(object):
	# constructor
	def __init__(self, discreteFunctions):
		# input class
		self.dF = discreteFunctions
		self.dFfullyInitialized = not (self.dF.a0Discrete == None
									   or self.dF.SigmaDiscrete == None
									   or self.dF.sigma0Discrete == None
									   or self.dF.Omega == None
									   or self.dF.kappa == None)
		# the W composition in powers of omega
		self.W0 = None
		self.W1 = None
		self.W2 = None
		self.W3 = None
		self.W4 = None
		self.W5 = None
		self.W = None

	# ===============
	# init functions
	# ===============

	def initW0(self):
		if self.W0 == None and self.dFfullyInitialized:
			self.W0 = self.__calcW0()
		else:
			raise BaseException("initW0(): already or DiscreteFunctions passed not fully initialized")

	def initW1(self):
		if self.W1 == None and self.dFfullyInitialized:
			self.W1 = self.__calcW1()
		else:
			raise BaseException("initW1(): already or DiscreteFunctions passed not fully initialized")

	def initW2(self):
		if self.W2 == None and self.dFfullyInitialized:
			self.W2 = self.__calcW2()
		else:
			raise BaseException("initW2(): already or DiscreteFunctions passed not fully initialized")

	def initW3(self):
		if self.W3 == None and self.dFfullyInitialized:
			self.W3 = self.__calcW3()
		else:
			raise BaseException("initW3(): already or DiscreteFunctions passed not fully initialized")

	def initW4(self):
		if self.W4 == None and self.dFfullyInitialized:
			self.W4 = self.__calcW4()
		else:
			raise BaseException("initW4(): already or DiscreteFunctions passed not fully initialized")

	def initW5(self):
		if self.W5 == None and self.dFfullyInitialized:
			self.W5 = self.__calcW5()
		else:
			raise BaseException("initW5(): already or DiscreteFunctions passed not fully initialized")

	def initW(self):
		if self.W == None \
				and not (self.W0 == None
						 or self.W1 == None
						 or self.W2 == None
						 or self.W3 == None
						 or self.W4 == None
						 or self.W5 == None)\
				and self.dFfullyInitialized:
			self.W = self.__calcW()
		else:
			raise BaseException("initW(): already, W# not or DiscreteFunctions passed not fully initialized")

	# ======================
	# calculating functions
	# ======================

	def __calcW0(self):
		pass

	def __calcW1(self):
		pass

	def __calcW2(self):
		pass

	def __calcW3(self):
		pass

	def __calcW4(self):
		pass

	def __calcW5(self):
		pass

	def __calcW(self):
		return self.W0 + self.W1 + self.W2 + self.W3 + self.W4 + self.W5

	# ==============================
	# component functions:
	# see calculations for reason
	# ==============================

	# todo: evaluate the kronecker delta by using if -> might be faster!

	def delta(self, i, j):
		res = None
		if i == j: res = 1
		else: res = 0
		return res

	def rValue(self, i):
		rVals = self.dF.dC.radialCells.rValues
		return rVals[i]

	def d1Value(self, i, j):
		return self.dF.linAlg.logDerivationMatrixOne[i][j]

	def d2Value(self, i, j):
		return self.dF.linAlg.logDerivationMatrixTwo[i][j]

	def jValue(self, i, j):
		return self.dF.linAlg.JMatrix[i][j]

	def einsum(self, function1, function2, i, k):
		res = 0.0
		for j in np.arange(self.dF.iP.number):
			res += function1(i, j) * function2(j, k)
		return res

	def fA0(self, i, k):
		Sigma = self.dF.SigmaDiscrete
		p, q = self.dF.iP.p, self.dF.iP.q
		summand1 = self.rValue(i) * self.einsum(self.d1Value, self.jValue, i, k)
		summand2 = self.d1Value(i, k)/(Sigma[i])
		summand3 = - q*self.delta(i, k)/(Sigma[i])
		summand4 = self.rValue(i)*(1-p)*self.einsum(self.delta, self.jValue, i, k)
		return summand1 + summand2 + summand3 + summand4

	def fA2(self,i, k):
		# constants
		m = self.dF.iP.m
		# todo: change to scipy constant
		G  = 6.67384*10**(-11) # in SI
		mStar = self.dF.iP.mStar
		mDisk = self.dF.iP.mDisk
		return self.delta(1, m) * self.rValue(i)**4/2*(self.jValue(i, k))/(G*(mStar+mDisk))

	def fB0(self, i, k):
		Sigma = self.dF.SigmaDiscrete
		summand1 = self.rValue(i)**2*self.einsum(self.delta, self.jValue, i, k)
		summand2 = self.rValue(i)*self.delta(i, k)/(Sigma[i])
		return summand1 + summand2

	def fB2(self, i, k):
		# constants
		m = self.dF.iP.m
		# todo: change to scipy constant
		G  = 6.67384*10**(-11) # in SI
		mStar = self.dF.iP.mStar
		mDisk = self.dF.iP.mDisk
		return self.delta(1, m) * self.rValue(i)**5/2*(self.jValue(i, k))/(G*(mStar+mDisk))

	def x0(self):
		pass

	def x1(self):
		pass

	def x2(self):
		pass

	def x3(self):
		pass

	def y0(self):
		pass

	def y1(self):
		pass

	def y2(self):
		pass

	def y3(self):
		pass

	def f1(self, i, k):
		Sigma = self.dF.SigmaDiscrete
		p, q = self.dF.iP.p, self.dF.iP.q
		summand1 = self.einsum(self.d2Value, self.jValue, i, k)
		summand2 = 2*(1-p)*self.einsum(self.d1Value, self.jValue, i, k)
		summand3 = - self.einsum(self.d1Value, self.jValue, i, k)
		summand4 = p*(p-1)*self.einsum(self.delta, self.jValue, i, k)
		return summand1 + summand2 + summand3 + summand4

	def f2(self):
		pass

	def b0(self):
		pass

	def b1(self):
		pass

	def b2(self):
		pass

	def b3(self):
		pass

	def c0(self):
		pass

	def c1(self):
		pass

	def c2(self):
		pass

	def c3(self):
		pass

	def c4(self):
		pass

	def c5(self):
		pass

"""

# =================== TESTING AREA ===================

import Params as P

iP = P.InitialParameters(number=5000)
dC = P.DerivedConstants(iP)

la = LinearAlgebraFunctions(iP, dC)
timeStart1 = t.clock()
la.initLDMOne()
la.initLDMTwo()
timeEnd1 = t.clock()
print "done1: took me [s]:", str(timeEnd1-timeStart1)
timeStart2 = t.clock()
la.initJMatrix()
timeEnd2 = t.clock()
print "done2: took me [s]:", str(timeEnd2-timeStart2)

"""
