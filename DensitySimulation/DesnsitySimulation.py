from __future__ import division, print_function

import time
import matplotlib
import numpy as np

import Params,Functions,LinAl

def densitySimulation(number=100, rStar=10 ** 11, rDisk=10 ** 15, p=1, q=5 / 8, mStar=10 ** 33, mDisk=10 ** 33, a0Star=10**9, m=1, eta=0.1):
	"""
	:param number: Number of cells to use. More means better resolution but slower run. Note cell size is on logarithmic scale.
	:param rStar: Radius of star in cm
	:param rDisk: Radius of disk in cm
	:param p: Power law index of density. Between 1/2 and 2. Probably better modern estimates available
	:param q: Power law index of temperature/sound speed. Between 1/2 and 3/4
	:param mStar: Mass of disk in grams
	:param mDisk: Mass of star in grams
	:param a0Star: speed of sound at rStar in cm/s
	:param m: Disk wave number, or mode
	:param eta: Gravity softening parameter
	:return:
	"""

	initParams=Params.InitialParameters(number,rStar,rDisk,p,q,mStar,mDisk,a0Star,m,eta)
	derivedConstants=Params.DerivedConstants(initParams)
	analyticFunctions=Functions.AnalyticalFunctions(initParams,derivedConstants)
	discreteFunctions=Functions.DiscreteFunctions(analyticFunctions)
	discreteFunctions.init()
	wMatrix=LinAl.WMatrix(discreteFunctions)
	wMatrix.init()
	eigenSolver=LinAl.EigenvalueSolver(wMatrix)
	eigenSolver.initEigen()

	print(eigenSolver.eigenvalues)


A=time.clock()
densitySimulation()
B=time.clock()
print(B-A)