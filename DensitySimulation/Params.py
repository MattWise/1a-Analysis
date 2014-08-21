from __future__ import division

import math as m
import numpy as np
from scipy import constants

""" Representation of a set of initial conditions.
"""
class InitialParameters(object):
	def __init__(self, number=500, rStar=10 ** 11, rDisk=10 ** 15, p=1, q=5 / 8, mStar=10 ** 33, mDisk=10 ** 33,
				 a0Star=10**6, m=1, eta=0.1):
		# do the assignments
		self.number = number
		self.rStar = rStar
		self.rDisk = rDisk
		self.p = p
		self.q = q
		self.mStar = mStar
		self.mDisk = mDisk
		self.a0Star = a0Star
		self.m = m
		self.eta = eta

	""" toString representation.
	"""
	def __str__(self):
		return "number: " + str(self.number) \
			   + "\n(rStar,rDisk): " + str((self.rStar, self.rDisk)) \
			   + "\n(mStar,mDisk): " + str((self.mStar, self.mDisk)) \
			   + "\n(p,q): " + str((self.p, self.q)) \
			   + "\na0Star: " + str(self.a0Star) \
			   + "\nm: " + str(self.m) \
			   + "\neta: " + str(self.eta)


""" Class for calculating all the derived constants out there. Based on the initial constants.
"""
class DerivedConstants(object):
	def __init__(self, initParams):
		# do the object creation.
		self.initialParameters = initParams
		self.n = 4 - (2 / initParams.q)
		self.sigmaStar = self.__sigStar(initParams)
		self.omegaStar = self.__omStar(initParams)
		self.qStar = self.__qStar(initParams)
		self.radialCells = self.__radCells(initParams)

	""" This class represents a set of radial grid points having given
	the number of cells and radius of star & radius of disk. Uses numpy arrays.
	"""
	class RadialCells(object):
		def __init__(self, number, rStar, rDisk):
			# simple constructor without setters, getters and stuff. thus, the user is
			# is responsible to use the class in a proper way.
			self.number = float(number)
			self.rStar = float(rStar)
			self.rDisk = float(rDisk)
			self.f = self.__calcF()
			self.rValues = None
			self.__calcRvalues()

		# ===============
		# inner methods:
		# ===============

		# definition see page 973, B1b
		def __calcF(self):
			return (self.rDisk / self.rStar) ** (1 / (self.number - 1))

		# definition see page 972, B1a
		def __calcRvalues(self):
			self.rValues = np.array([])
			for i in range(1, int(self.number)):
				value = self.f ** (i - 1) * self.rStar
				self.rValues = np.append(self.rValues, value)
			# here, we assert that r_N = rDisk for completness but
			# due to numerical errors, we have to fix that the following
			# way
			self.rValues = np.append(self.rValues, self.rDisk)

		# simple toString method for proper debug and information
		def __str__(self):
			return "number: " + str(self.number) + "\n" \
				   + "rStar: " + str(self.rStar) + "\n" \
				   + "rDisk: " + str(self.rDisk)

	""" toString representation.
	"""
	def __str__(self):
		return "n: " + str(self.n) + "\n" \
			   + "sigmaStar: " + str(self.sigmaStar) + "\n" \
			   + "omegaStar: " + str(self.omegaStar) + "\n" \
			   + "qStar: " + str(self.qStar)

	""" sigmaStar constant for analytic form of sigma0. definition page 964, 23a and 23b.
	"""
	def __sigStar(self, iP):
		return (2 - iP.p) * iP.mStar / \
			   (2 * m.pi * iP.rStar ** 2 * ((iP.rDisk / iP.rStar) ** (2 - iP.p) - 1))

	""" omegaStar constant for qStar constant. definition page 964, 24.
	"""
	def __omStar(self, iP):
		return m.sqrt(constants.G * iP.mStar / (iP.rStar ** 2))

	""" qStar constant. definition page 964, 24.
	"""
	def __qStar(self, iP):
		return (self.omegaStar * iP.a0Star) / (m.pi * constants.G * self.sigmaStar)

	""" creates a representation of the grid in respect to the initial parameter n.
	"""
	def __radCells(self, iP):
		return self.RadialCells(iP.number, iP.rStar, iP.rDisk)