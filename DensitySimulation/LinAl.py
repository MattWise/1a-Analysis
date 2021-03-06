from __future__ import division

import numpy as np
from scipy.constants import pi, G
import scipy.linalg as sciLa

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
		print "-> begin building J matrix"
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
			# todo: do it sideways for better overview
			percentage = round(100*(i/dimension), 0)
			if percentage%5==0 and not (percentage==printed):
				print "--> completed:", percentage, "%"
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
									   or self.dF.gSigma0Discrete == None
									   or self.dF.Omega == None
									   or self.dF.kappa == None)
		# access to important fields
		self.m = self.dF.iP.m
		self.p = self.dF.iP.p
		self.q = self.dF.iP.q
		self.mRatio = self.dF.iP.mRatio
		self.number = self.dF.iP.number
		# the W composition in powers of omega
		self.W0 = None
		self.W1 = None
		self.W2 = None
		self.W3 = None
		self.W4 = None
		self.W5 = None
		self.B14A = None
		self.B14C = None

	# ===============
	# init functions
	# ===============

	def init(self):
		# todo: parallel?
		self.initW0()
		self.initW1()
		self.initW2()
		self.initW3()
		self.initW4()
		self.initW5()
		self.initB14A()
		self.initB14C()

	def initW0(self):
		if self.W0 == None and self.dFfullyInitialized:
			print("initW0(): start")
			self.W0 = self.__calcW0()
			print("initW0(): end")
		else:
			raise BaseException("initW0(): already or DiscreteFunctions passed not fully initialized")

	def initW1(self):
		if self.W1 == None and self.dFfullyInitialized:
			print("initW1(): start")
			self.W1 = self.__calcW1()
			print("initW1(): end")
		else:
			raise BaseException("initW1(): already or DiscreteFunctions passed not fully initialized")

	def initW2(self):
		if self.W2 == None and self.dFfullyInitialized:
			print("initW2(): start")
			self.W2 = self.__calcW2()
			print("initW2(): end")
		else:
			raise BaseException("initW2(): already or DiscreteFunctions passed not fully initialized")

	def initW3(self):
		if self.W3 == None and self.dFfullyInitialized:
			print("initW3(): start")
			self.W3 = self.__calcW3()
			print("initW3(): end")
		else:
			raise BaseException("initW3(): already or DiscreteFunctions passed not fully initialized")

	def initW4(self):
		if self.W4 == None and self.dFfullyInitialized:
			print("initW4(): start")
			self.W4 = self.__calcW4()
			print("initW4(): end")
		else:
			raise BaseException("initW4(): already or DiscreteFunctions passed not fully initialized")

	def initW5(self):
		if self.W5 == None and self.dFfullyInitialized:
			print("initW5(): start")
			self.W5 = self.__calcW5()
			print("initW5(): end")
		else:
			raise BaseException("initW5(): already or DiscreteFunctions passed not fully initialized")

	def initB14A(self):
		if self.B14A == None \
				and not (self.W0 == None
						 or self.W1 == None
						 or self.W2 == None
						 or self.W3 == None
						 or self.W4 == None
						 or self.W5 == None)\
				and self.dFfullyInitialized:
			print("initB14A(): start")
			self.B14A = self.__calcB14A()
			print("initB14A(): end")
		else: raise BaseException("initB14A(): already, W# not or DiscreteFunctions passed not fully initialized")


	def initB14C(self):
		if self.B14C == None \
				and not (self.W0 == None
						 or self.W1 == None
						 or self.W2 == None
						 or self.W3 == None
						 or self.W4 == None
						 or self.W5 == None)\
				and self.dFfullyInitialized:
			print("initB14C(): start")
			self.B14C = self.__calcB14C()
			print("initB14C(): end")
		else: raise BaseException("initB14C(): already, W# not or DiscreteFunctions passed not fully initialized")

	# ======================
	# calculating functions
	# ======================

	def __calcW0(self):
		res = np.zeros((self.number, self.number), dtype=np.complex128)
		for i in np.arange(self.number):
			for k in np.arange(self.number):
				res[i][k] = self.calcW0_ik(i, k)
		return res

	def __calcW1(self):
		res = np.zeros((self.number, self.number), dtype=np.complex128)
		for i in np.arange(self.number):
			for k in np.arange(self.number):
				res[i][k] = self.calcW1_ik(i, k)
		return res

	def __calcW2(self):
		res = np.zeros((self.number, self.number), dtype=np.complex128)
		for i in np.arange(self.number):
			for k in np.arange(self.number):
				res[i][k] = self.calcW2_ik(i, k)
		return res

	def __calcW3(self):
		res = np.zeros((self.number, self.number), dtype=np.complex128)
		for i in np.arange(self.number):
			for k in np.arange(self.number):
				res[i][k] = self.calcW3_ik(i, k)
		return res

	def __calcW4(self):
		res = np.zeros((self.number, self.number), dtype=np.complex128)
		for i in np.arange(self.number):
			for k in np.arange(self.number):
				res[i][k] = self.calcW4_ik(i, k)
		return res

	def __calcW5(self):
		res = np.zeros((self.number, self.number), dtype=np.complex128)
		for i in np.arange(self.number):
			for k in np.arange(self.number):
				res[i][k] = self.calcW5_ik(i, k)
		return res

	def __calcB14A(self):
		zeros=np.zeros((self.number,self.number),dtype=np.complex128)
		I=np.identity(self.number,dtype=np.complex128)
		C1=np.concatenate((-1*self.W0,zeros,zeros,zeros,zeros),0)
		C2=np.concatenate((zeros,I,zeros,zeros,zeros),0)
		C3=np.concatenate((zeros,zeros,I,zeros,zeros),0)
		C4=np.concatenate((zeros,zeros,zeros,I,zeros),0)
		C5=np.concatenate((zeros,zeros,zeros,zeros,I),0)
		columns=(C1,C2,C3,C4,C5)
		array=np.concatenate(columns,1)
		return array

	def __calcB14C(self):
		zeros=np.zeros((self.number,self.number),dtype=np.complex128)
		I=np.identity(self.number,dtype=np.complex128)
		C1=np.concatenate((self.W1,self.W5,self.W4,self.W3,self.W2),0)
		C2=np.concatenate((zeros,zeros,I,zeros,zeros),0)
		C3=np.concatenate((zeros,zeros,zeros,I,zeros),0)
		C4=np.concatenate((zeros,zeros,zeros,zeros,I),0)
		C5=np.concatenate((I,zeros,zeros,zeros,zeros),0)
		columns=(C1,C2,C3,C4,C5)
		array=np.concatenate(columns,1)
		return array


	# ==============================
	# value functions
	# ==============================

	def delta(self, i, j):
		return int(i==j)

	def rValue(self, i):
		return self.dF.dC.radialCells.rValues[i]

	def d1Value(self, i, j):
		return self.dF.linAlg.logDerivationMatrixOne[i][j]

	def d2Value(self, i, j):
		return self.dF.linAlg.logDerivationMatrixTwo[i][j]

	def jValue(self, i, j):
		return self.dF.linAlg.JMatrix[i][j]

	def SigmaValue(self, i):
		return self.dF.SigmaDiscrete[i]

	def gSigma0Value(self, i):
		return self.dF.gSigma0Discrete[i]

	def DgSigma0Value(self,i):
		return self.dF.DgSigma0Discrete[i]

	def kappaValue(self, i):
		return self.dF.kappa[i]

	def DkappaValue(self,i):
		return self.dF.Dkappa[1]

	def OmegaValue(self, i):
		return self.dF.Omega[i]

	def DOmegaValue(self,i):
		return self.dF.DOmega[i]

	def DDOmegaValue(self,i):
		return self.dF.DDOmega(i)

	# einstein sum
	def einsum(self, function1, function2, i, k):
		res = 0.0
		for j in np.arange(self.dF.iP.number):
			res += function1(i, j) * function2(j, k)
		return res

	# ======================================
	# component functions:
	# see calculations for further detail
	# ======================================
	# todo: evaluate the kronecker delta by using if -> might be faster but worse to read.

	def fA0(self, i, k):
		summand1 = self.rValue(i) * self.einsum(self.d1Value, self.jValue, i, k)
		summand2 = self.d1Value(i, k)/(self.gSigma0Value(i))
		summand3 = - self.q*self.delta(i, k)/(self.gSigma0Value(i))
		summand4 = self.rValue(i)*(1-self.p)*self.einsum(self.delta, self.jValue, i, k)
		return summand1 + summand2 + summand3 + summand4

	def fA2(self,i, k):
		return self.delta(1, self.m) * self.rValue(i)**4/2\
			   *(self.jValue(i, k))/(1+self.mRatio)

	def fB0(self, i, k):
		summand1 = self.rValue(i)**2*self.einsum(self.delta, self.jValue, i, k)
		summand2 = self.rValue(i)*self.delta(i, k)/(self.SigmaValue(i))
		return summand1 + summand2

	def fB2(self, i, k):
		"""
		if (i in range(6,9)) and (k in range(0,10)):
			print "rValue=", self.rValue(i), "---> rValue**5=", self.rValue(i)**5
			print "jValue", self.jValue(i, k)
			print "mRatio", self.mRatio
		"""
		return self.delta(1, self.m) * self.rValue(i)**5/2\
			   *(self.jValue(i, k))/(1+self.mRatio)

	def x0(self, i, k):
		denominator=self.rValue(i)*self.kappaValue(i)**3*self.gSigma0Value(i)
		summand1 = self.DgSigma0Value(i)*self.rValue(i)*(self.m**2*self.OmegaValue(i)**2-self.kappaValue(i)**2)
		summand2 = self.gSigma0Value(i)*\
					(\
						-self.kappaValue(i)**2\
						+2*self.kappaValue(i)*self.DkappaValue(i)*self.rValue(i)\
						+self.m**2*self.OmegaValue(i)*(self.OmegaValue(i)-2*self.DOmegaValue(i)*self.rValue(i))
					)
		return self.m*self.OmegaValue(i)*(\
			summand1 + summand2\
			)/denominator

	def x1(self, i, k):
		denominator=self.rValue(i)*self.kappaValue(i)**3*self.gSigma0Value(i)
		summand1 = self.DgSigma0Value(i)*(self.kappaValue(i)**2-3*self.m**2*self.OmegaValue(i)**2)*self.rValue(i)
		summand2 = self.gSigma0Value(i)*(\
			self.kappaValue(i)**2\
		    -2*self.kappaValue(i)*self.DkappaValue(i)*self.rValue(i)
		    +self.m**2*self.OmegaValue(i)*(-3*self.OmegaValue(i)+4*self.DOmegaValue(i)*self.rValue(i))
			)
		return (summand1 + summand2)/denominator

	def x2(self, i, k):
		denominator=self.rValue(i)*self.kappaValue(i)**3*self.gSigma0Value(i)
		summand1 = 3*self.DgSigma0Value(i)*self.OmegaValue(i)*self.rValue(i)
		summand2 = 3*self.OmegaValue(i)*self.gSigma0Value(i)
		summand3 = -2*self.DOmegaValue(i)*self.rValue(i)*self.gSigma0Value(i)
		return self.m*(summand1 + summand2 + summand3)/denominator

	def x3(self, i, k):
		numerator=self.gSigma0Value(i)+self.rValue(i)*self.DgSigma0Value(i)
		denominator=self.rValue(i)*self.kappaValue(i)**3*self.gSigma0Value(i)
		return -numerator/denominator

	def y0(self, i, k):
		summand1 = 2*self.DgSigma0Value(i)*self.OmegaValue(i)*self.rValue(i)\
		           *(self.m**2*self.OmegaValue(i)**2-self.kappaValue(i)**2)
		summand2 = self.gSigma0Value(i)*(\
			4*self.DkappaValue(i)*self.kappaValue(i)*self.OmegaValue(i)*self.rValue(i)\
			+ self.kappaValue(i)**2*(self.m**2*self.OmegaValue(i)-2*self.DOmegaValue(i)*self.rValue(i))\
			- self.m**2*self.OmegaValue(i)**2*(self.m**2*self.OmegaValue(i)+2*self.rValue(i)*self.DOmegaValue(i))\
			)
		denominator=self.rValue(i)**2*self.kappaValue(i)**3*self.gSigma0Value(i)
		return self.m*(summand1 + summand2)/denominator

	def y1(self, i, k):
		term1 = self.m**2
		term2part1 = -self.kappaValue(i)**2*self.gSigma0Value(i)
		term2part2 = 3*self.m**2*self.gSigma0Value(i)*self.OmegaValue(i)**2
		term2part3 = -4*self.rValue(i)*self.OmegaValue(i)**2*self.DgSigma0Value(i)
		denominator = self.rValue(i)**2*self.kappaValue(i)**3*self.gSigma0Value(i)
		return term1*(term2part1+term2part2+term2part3)/denominator

	def y2(self, i, k):
		term1 = 2*self.DgSigma0Value(i)*self.OmegaValue(i)*self.rValue(i)
		term2 = -3*self.m**2*self.OmegaValue(i)*self.gSigma0Value(i)
		term3 = 2*self.DOmegaValue(i)*self.rValue(i)*self.gSigma0Value(i)
		denominator=self.rValue(i)**2*self.kappaValue(i)**3*self.gSigma0Value(i)
		return self.m*(term1+term2+term3)/denominator

	def y3(self, i, k):
		numerator=self.m**2
		denominator=self.rValue(i)**2*self.kappaValue(i)**3
		return numerator/denominator

	def f1(self, i, k):
		summand1 = self.einsum(self.d2Value, self.jValue, i, k)
		summand2 = 2*(1-self.p)*self.einsum(self.d1Value, self.jValue, i, k)
		summand3 = - self.einsum(self.d1Value, self.jValue, i, k)
		summand4 = self.p*(self.p-1)*self.einsum(self.delta, self.jValue, i, k)
		summand5 = self.d2Value(i, k)/(self.rValue(i)*self.gSigma0Value(i))
		summand6 = - 2*self.q*self.d1Value(i, k)/(self.rValue(i)*self.gSigma0Value(i))
		summand7 = - self.d1Value(i, k)/(self.rValue(i)*self.gSigma0Value(i))
		summand8 = self.q*(self.q+1)*self.delta(i, k)/(self.rValue(i)*self.gSigma0Value(i))
		return summand1 \
			   + summand2 \
			   + summand3 \
			   + summand4 \
			   + summand5 \
			   + summand6 \
			   + summand7 \
			   + summand8

	def f2(self, i, k):
		return - (self.kappaValue(i)**2*self.rValue(i)*self.delta(i, k))/(2*pi*self.gSigma0Value(i))

	def b0(self, i, k):
		summand1 = - self.m * self.OmegaValue(i)*self.kappaValue(i)**2
		summand2 = self.m**3 * self.OmegaValue(i)**3
		return summand1 + summand2

	def b1(self, i, k):
		summand1 = self.kappaValue(i)**2
		summand2 = - 3 * self.m**2*self.OmegaValue(i)**2
		return summand1 + summand2

	def b2(self, i, k):
		return 3 * self.m * self.OmegaValue(i)

	def b3(self, i, k):
		return -1

	def c0(self, i, k):
		summand1 = 2 * self.m**3 * self.kappaValue(i)**2 * self.OmegaValue(i)**3
		summand2 = - self.m**5 * self.OmegaValue(i)**5
		summand3 = - self.m * self.kappaValue(i)**4 * self.OmegaValue(i)
		return summand1 + summand2 + summand3

	def c1(self, i, k):
		summand1 = self.kappaValue(i)**4
		summand2 = - 6 * self.m**2 * self.kappaValue(i)**2 * self.OmegaValue(i)**2
		summand3 = 5 * self.m**4 * self.OmegaValue(i)**4
		return summand1 + summand2 + summand3

	def c2(self, i, k):
		summand1 = 6 * self.m * self.OmegaValue(i) * self.kappaValue(i)**2
		summand2 = 10 * self.m**3 * self.OmegaValue(i)**3
		return summand1 + summand2

	def c3(self, i, k):
		summand1 = 2 * self.kappaValue(i)**2
		summand2 = - 2 * 5 * self.m**2 * self.OmegaValue(i)**2
		return summand1 + summand2

	def c4(self, i, k):
		return - 5 * self.m * self.OmegaValue(i)

	def c5(self, i, k):
		return 1


	# ===========================================
	# functions for calculating W^(x)_ik element
	# ===========================================

	# todo: make parallel :-O uses less memory and only one processor at 100%

	def calcW0_ik(self, i, k):
		if i==0:
			summand1 = -self.m*self.OmegaValue(i)*self.einsum(self.d1Value, self.jValue, i, k)
			summand2 = -self.m*self.OmegaValue(i)*(1-self.p)*self.einsum(self.delta, self.jValue, i, k)
			summand3 = -2*self.m*self.OmegaValue(i)*self.einsum(self.delta, self.jValue, i, k)
			summand4 = -self.m*self.OmegaValue(i)*self.d1Value(i, k)/(self.SigmaValue(i)*self.rValue(i))
			summand5 = self.m*self.OmegaValue(i)*self.q*self.delta(i, k)/(self.SigmaValue(i)*self.rValue(i))
			summand6 = -2*self.m*self.OmegaValue(i)*self.delta(i, k)/(self.SigmaValue(i)*self.rValue(i))
			return summand1 + summand2 + summand3+ summand4 + summand5 + summand6
		elif i==self.number-1:
			summand1 = -self.m*self.OmegaValue(i)*self.einsum(self.d1Value, self.jValue, i, k)
			summand2 = -self.m*self.OmegaValue(i)*(1-self.p)*self.einsum(self.delta, self.jValue, i, k)
			summand3 = -2*self.m*self.OmegaValue(i)*self.einsum(self.delta, self.jValue, i, k)
			summand4 = -self.m*self.OmegaValue(i)*self.d1Value(i, k)/(self.SigmaValue(i)*self.rValue(i))
			summand5 = self.m*self.OmegaValue(i)*self.q*self.delta(i, k)/(self.SigmaValue(i)*self.rValue(i))
			summand6 = -2*self.m*self.OmegaValue(i)*self.delta(i, k)/(self.SigmaValue(i)*self.rValue(i))
			summand7 = - self.rValue(i)*self.delta(i, k)*self.m\
					   *self.OmegaValue(i)/(2*pi*self.gSigma0Value(i)*self.p)\
					   *(self.kappaValue(i)**2-self.m**2*self.OmegaValue(i)**2)
			return summand1 + summand2 + summand3 \
					+ summand4 + summand5 + summand6 + summand7
		else:
			summand1 = self.fA0(i, k) * (self.x0(i, k) + self.x2(i, k))
			summand2 = self.fB0(i, k) * (self.y0(i, k) + self.y2(i, k))
			summand3 = (self.c0(i, k) * self.f2(i, k))/(self.kappaValue(i)**5)
			summand4 = (self.b0(i, k) * self.f1(i, k))/(self.kappaValue(i)**3)
			return summand1 + summand2 + summand3 + summand4

	def calcW1_ik(self, i, k):
		if i==0:
			summand1 = self.einsum(self.d1Value, self.jValue, i, k)
			summand2 = (1-self.p)*self.einsum(self.delta, self.jValue, i, k)
			summand3 = - self.delta(i, k) * self.q/(self.SigmaValue(i)*self.rValue(i))
			summand4 = self.d1Value(i, k)/(self.SigmaValue(i)*self.rValue(i))
			return summand1 + summand2 + summand3 + summand4
		elif i==self.number-1:
			summand1 = self.einsum(self.d1Value, self.jValue, i, k)
			summand2 = (1-self.p)*self.einsum(self.delta, self.jValue, i, k)
			summand3 = self.d1Value(i, k)/(self.SigmaValue(i)*self.rValue(i))
			summand4 = - self.q * self.delta(i, k)/(self.SigmaValue(i)*self.rValue(i))
			summand5 = self.rValue(i)*self.delta(i, k)\
					   *(self.kappaValue(i)**2-self.m**2*self.OmegaValue(i)**2)/(2*pi*self.p*self.gSigma0Value(i))
			summand6 = - self.rValue(i)*self.delta(i, k) \
					   * self.m**2*self.OmegaValue(i)**2/(2*pi*self.gSigma0Value(i)*self.p)
			return summand1 + summand2 + summand3 + summand4 + summand5 + summand6
		else:
			summand1 = self.fA0(i, k) * self.x1(i, k)
			summand2 = self.fB0(i, k) * self.y1(i, k)
			summand3 = (self.c1(i, k) * self.f2(i, k))/(self.kappaValue(i)**5)
			summand4 = (self.b1(i, k) * self.f1(i, k))/(self.kappaValue(i)**3)
			return summand1 + summand2 + summand3 + summand4

	def calcW2_ik(self, i, k):
		if i==0:
			return - self.delta(1, self.m) *3/2*self.OmegaValue(i)*self.rValue(i)**3\
				   *self.jValue(i, k)/(1+self.mRatio)
		elif i==self.number-1:
			summand1 = -self.delta(1, self.m)*3*self.OmegaValue(i)*self.rValue(i)**3\
					   *self.jValue(i, k)/(2*(1+self.mRatio))
			summand2 = - self.m * self.OmegaValue(i) * self.rValue(i) * (-1) \
					   * self.delta(i, k)/(2*pi*self.p*self.gSigma0Value(i))
			summand3 = self.rValue(i)*self.delta(i, k)*2*self.m*self.OmegaValue(i)\
					   /(2*pi*self.p*self.gSigma0Value(i))
			return summand1 + summand2 + summand3
		else:
			summand1 = self.fA0(i, k)
			summand2 = self.fB0(i, k)
			summand3 = self.fA2(i, k) * (self.x0(i, k) + self.x2(i, k))
			summand4 = self.fB2(i, k) * (self.y0(i, k) + self.y2(i, k))
			summand5 = (self.c2(i, k) * self.f2(i, k))/(self.kappaValue(i)**5)
			summand6 = (self.b2(i, k) * self.f1(i, k))/(self.kappaValue(i)**3)
			return summand1 + summand2 + summand3 + summand4 + summand5 + summand6

	def calcW3_ik(self, i, k):
		if i==0:
			return self.delta(1, self.m)*self.rValue(i)**3*self.jValue(i, k)/(2*(1+self.mRatio))
		elif i==self.number-1:
			summand1 = (self.delta(1, self.m)*self.rValue(i)**3*self.jValue(i, k))/(2*(1+self.mRatio))
			summand2 = (self.rValue(i)*self.delta(i, k)*(-1))/(2*pi*self.gSigma0Value(i)*self.p)
			return summand1 + summand2
		else:
			summand1 = self.fA2(i, k) * self.x1(i, k)
			summand2 = self.fA0(i, k) * self.x3(i, k)
			summand3 = self.fB2(i, k) * self.y1(i, k)
			summand4 = self.fB0(i, k) * self.y3(i, k)
			summand5 = (self.c3(i, k) * self.f2(i, k))/(self.kappaValue(i)**3)
			summand6 = (self.b3(i, k) * self.f1(i, k))/(self.kappaValue(i)**3)
			return summand1 + summand2 + summand3 + summand4 + summand5 + summand6

	def calcW4_ik(self, i, k):
		if i==0:
			return 0
		elif i==self.number-1:
			return 0
		else:
			summand1 = self.fA2(i, k)
			summand2 = self.fB2(i, k)
			summand3 = (self.c4(i, k) * self.f2(i, k))/(self.kappaValue(i)**5)
			return summand1 + summand2 + summand3

	def calcW5_ik(self, i, k):
		if i==0:
			return 0
		elif i==self.number-1:
			return 0
		else:
			summand1 = self.fA2(i, k) * self.x3(i, k)
			summand2 = self.fB2(i, k) * self.y3(i, k)
			summand3 = (self.c5(i, k) * self.f2(i, k))/(self.kappaValue(i)**5)
			"""
			if i==8 and (k in range(0,10)):
				print "======"
				print (i,k), "summand1", summand1
				print (i,k), "summand2", summand2
				print "--summand2--->"
				print "------fB2---->", self.fB2(i,k)
				print "------y3----->", self.y3(i,k)
				print (i,k), "summand3", summand3
			"""
			return summand1 + summand2 + summand3


""" Given an instance of a fully initialized WMatrix object, this class uses the matrices B14A and B14C. The result is
the wanted S vector (surface density purturbance).
"""
class EigenvalueSolver(object):
	# constructor
	def __init__(self, wMatrix):
		# input class
		self.wM = wMatrix
		self.eigenvalues = None
		self.eigenvectors = None
		self.S = None

	def initEigen(self):
		if self.wM.dFfullyInitialized \
				and not (self.wM.B14A == None or self.wM.B14C == None):
			self.eigenvalues, self.eigenvectors = self.__calcEigen()
		else:
			raise BaseException("given WMatrix object not fully or its B14A, B14C matrices not initialized")
	
	def __calcEigen(self):
		return sciLa.eig(self.wM.B14A, self.wM.B14C)