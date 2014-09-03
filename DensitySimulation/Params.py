from __future__ import division

import math as m
import numpy as np

""" Representation of a set of initial conditions.
"""
class InitialParameters(object):
	def __init__(self, number=500, rRatio=10**4, p=1, q=5 / 8, mRatio=1,
				 qStar=.3, m=1, eta=0.1):
		# do the assignments
		self.number = number
		self.rRatio = rRatio # equal to rDisk/rStar. Since rStar=1 because of natural units, it is equivalent to rDisk
		self.p = p
		self.q = q
		self.mRatio = mRatio
		self.qStar = qStar
		self.m = m
		self.eta = eta

	""" toString representation.
	"""
	def __str__(self):
		return "number: " + str(self.number) \
			   + "\nrRatio: " + str(self.rRatio) \
			   + "\nmRatio: " + str(self.mRatio) \
			   + "\n(p,q): " + str((self.p, self.q)) \
			   + "\nqStar: " + str(self.qStar) \
			   + "\nm: " + str(self.m) \
			   + "\neta: " + str(self.eta)


""" Class for calculating all the derived constants out there. Based on the initial constants.
"""
class DerivedConstants(object):
	def __init__(self, initParams):
		# do the object creation.
		self.initialParameters = initParams
		self.n = 4 - (2 / initParams.q)
		self.gSigma0Star = self.__gSig0Star(initParams)
		self.a0Star = self.__a0Star(initParams)
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
			return "[number: " + str(self.number) + " -- " \
				   + "rStar: " + str(self.rStar) + " -- " \
				   + "rDisk: " + str(self.rDisk) + "]"

	""" toString representation.
	"""
	def __str__(self):
		return "n: " + str(self.n) + "\n" \
			   + "gSigma0Star: " + str(self.gSigma0Star) + "\n" \
			   + "a0Star: " + str(self.a0Star) + "\n" \
			   + "radialCells: " + str(self.radialCells)

	""" sigmaStar constant for analytic form of sigma0. definition page 964, 23a and 23b.
	"""
	def __gSig0Star(self, iP):
		return (2 - iP.p) * (self.initialParameters.mRatio) / \
			   (2 * m.pi * ((iP.rRatio) ** (2 - iP.p) - 1))

	""" aStar constant calculated from the other constants. derived, used 964, 24 and 964, 23b.
	"""
	def __a0Star(self, iP):
		return iP.qStar*(2-iP.p)*iP.mRatio/(2*(iP.rRatio**(2-iP.p)-1))

	""" creates a representation of the grid in respect to the initial parameter n.
	"""
	def __radCells(self, iP):
		return self.RadialCells(iP.number, 1, iP.rRatio)