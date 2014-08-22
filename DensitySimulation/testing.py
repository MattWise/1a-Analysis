from __future__ import division,print_function

import Functions, LinAl, Params
import numpy as np
import time as t

def timing(function):
	timeStart = t.clock()
	function()
	timeEnd = t.clock()
	print("took me:", (timeEnd-timeStart), "s")

# todo: number=10 -> complexWarning, unclear why.
iP = Params.InitialParameters(number=100)
dC = Params.DerivedConstants(iP)

anaFcts = Functions.AnalyticalFunctions(iP, dC)

dsctFcts = Functions.DiscreteFunctions(anaFcts)
dsctFcts.initLDMOne()
dsctFcts.initLMDTwo()
dsctFcts.initJMatrix()

# todo: complex warning in _quadpack._qagse !
dsctFcts.initOmega()

dsctFcts.initDOmega()

dsctFcts.initDDOmega()

dsctFcts.initKappa()

# todo: complex warning in array2[index]=[...] !
dsctFcts.initDkappa()

dsctFcts.initA0Discrete()

dsctFcts.initSigmaDiscrete()

dsctFcts.initSigma0Discrete()

dsctFcts.initDSigma0Discrete()


""" debug println's

print("Omega:", dsctFcts.Omega)
print("DOmega:", dsctFcts.DOmega)
print("DDOmega:", dsctFcts.DDOmega)
print("kappa:", dsctFcts.kappa)
print("Dkappa:", dsctFcts.Dkappa)
print("a0:", dsctFcts.a0Discrete)
print("Sigma:", dsctFcts.SigmaDiscrete)
print("sigma0:", dsctFcts.sigma0Discrete)
print("Dsigma0:", dsctFcts.Dsigma0Discrete)
print("fully initialized:", wM.dFfullyInitialized)


"""


wM = LinAl.WMatrix(dsctFcts)

def init():
	wM.initW0()
	wM.initW1()
	wM.initW2()
	wM.initW3()
	wM.initW4()
	wM.initW5()

timing(init)


print("W0:", wM.W0)
print("W1:", wM.W1)
print("W2:", wM.W2)
print("W3:", wM.W3)
print("W4:", wM.W4)
print("W5:", wM.W5)