from __future__ import division,print_function

import numpy as np
import scipy as sp
from scipy import constants
import math as m
import params

class RadialCells:

	""" This class represents a set of radial grid points having given
	the intended number of cells for radius of star and radius of disk.
	"""
	def __init__(self, number, rStar, rDisk):
		self.number = float(number)
		self.rStar = float(rStar)
		self.rDisk = float(rDisk)
		self.f = self.__calcF()
		self.rValues = []
		self.__calcRvalues()

	# not to be called from the outside
	def __calcF(self):
		return (self.rDisk/self.rStar)**(1/(self.number-1))

	# not to be called from the outside;
	def __calcRvalues(self):
		for i in range(1, int(self.number)):
			self.rValues.append(self.f**(i-1) * self.rStar)
		self.rValues.append(self.rDisk)

	def __str__(self):
		return "number: " + str(self.number) + "\n" \
				+ "rStar: " + str(self.rStar) + "\n" \
				+ "rDisk: " + str(self.rDisk)

class LinearAlgebra:

	def __init__(self,parameters,derived_constants):
		self.parameters=parameters
		self.derived_constants=derived_constants
		self.generate_matrices()


	def differential_array_maker(self,order,inner_boundary=False,outer_boundary=False,size=50):

		#creates a differentiation matrix for numerical differentiation. "Order=1" means first derivative, "order=2" means second.

		def Diagnols(array,a,b,value):
			#goes down diagonal form specified coordinates, setting each to value.
			if a>b:
				times = size-a-1
			else:
				times=size-b
			for row in range(times):
				array[a,b]=value
				a+=1
				b+=1
			return array

		def SetBoundary(array,size,inner_boundary,outer_boundary):
			#sets the matrix boundary conditions. Defaults are from ARS paper.

			outer_boundary.reverse()

			for index,value in enumerate(inner_boundary):
				array[0,index]=value
			for index,value in enumerate(outer_boundary):
				array[size-1,size-1-index]=value
			return array

		array = np.zeros([size,size])
		if order == 1:
			array = Diagnols(array,1,0,-1)
			array = Diagnols(array,1,2,1)

			#defualt boundary conditions
			if not inner_boundary:
				inner_boundary=[-3,4,-1]
			if not outer_boundary:
				outer_boundary=[1,-4,3]

		if order == 2:
			array = Diagnols(array,1,0,1)
			array = Diagnols(array,1,1,-2)
			array = Diagnols(array,1,2,1)

			#defualt boundary conditions
			if not inner_boundary:
				inner_boundary = [2,-5,4,-1]
			if not outer_boundary:
				outer_boundary=[-14,-5,2]

		array=SetBoundary(array,size,inner_boundary,outer_boundary)

		return array

	def generate_matrices(self):
		self.D1=(1/(2*np.log10(self.derived_constants.radial_cells.f)))*self.differential_array_maker(1,size=self.parameters.N)
		self.D2=(1/(np.log10(self.derived_constants.radial_cells.f))**2)*self.differential_array_maker(2,size=self.parameters.N)
		self.I=np.linalg.inv(self.D1)

	def first_derivative(self,vector,r):
		return (1/r)*np.dot(self.D1,vector)

	def second_derivative(self,vector,r):
		return (1/r)**2*(np.dot(self.D2,vector)-np.dot(self.D1,vector))

	def integrate(self,vector,r):
		#Not right! Missing some factor, like the (1/r) in the first derivative.
		return np.dot(self.I,vector)

def Discretize(func,derived_parameters):
	#takes continuous function of r as input and outputs a vector with component j equal to value at jth cell
	R_List=derived_parameters.radial_cells.rValues
	Value_List=[]
	for value in R_List:
		Value_List.append(func(value))
	return np.array(Value_List)