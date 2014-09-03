from __future__ import division, print_function

import time
import matplotlib
import numpy as np

import Params,Functions,LinAl,Correctness

def densitySimulation(number=500, rRatio=10**4, p=1, q=5 / 8, mRatio=1, qStar=.3, m=1, eta=0.1):
	"""
	:param number: Number of cells to use. More means better resolution but slower run. Note cell size is on logarithmic scale.
	:param rRatio: Radius of the disk divided by the radius of the star
	:param p: Power law index of density. Between 1/2 and 2. Probably better modern estimates available
	:param q: Power law index of temperature/sound speed. Between 1/2 and 3/4
	:param mRatio: Mass of the disk divided by the mass of the star
	:param qStar: Toomre Q at the radius of the star
	:param m: Disk wave number, or mode
	:param eta: Gravity softening parameter
	:return:
	"""

	initParams=Params.InitialParameters(number,rRatio,p,q,mRatio,qStar,m,eta)
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
densitySimulation(number=20)
B=time.clock()
print(B-A)