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
		# integration matrix for J
		# todo: unclear if this is the way to go
		self.I = None

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

	# todo: not sure if something to keep
	def integrate(self,vector):
		#Not right! Missing some factor, like the (1/r) in the first derivative.
		return np.dot(self.I,vector)

	# todo: I'd like to have a matrix called J for convenient access.


""" This class represents the huge matrix generated in order to solve the omega eigenvalue and
the S eigenvector.
Due to heavy computing, needs to be initialized as well.
"""
class WMatrix(object):
	pass