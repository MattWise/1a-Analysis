from __future__ import division,print_function

import Functions, LinAl, Params
import numpy as np
import time
import multiprocessing

np.set_printoptions(threshold=np.nan,linewidth=2000)

"""
iP = Params.InitialParameters(number=100)
print(str(iP))
dC = Params.DerivedConstants(iP)
print(str(dC))

anaFcts = Functions.AnalyticalFunctions(iP, dC)




dsctFcts = Functions.DiscreteFunctions(anaFcts)

dsctFcts.initOmega()
print("Omega:", dsctFcts.Omega)

dsctFcts.initDOmega()
print("DOmega:", dsctFcts.DOmega)

dsctFcts.initDDOmega()
print("DDOmega:", dsctFcts.DDOmega)

dsctFcts.initKappa()
print("kappa:", dsctFcts.kappa)

dsctFcts.initDkappa()
print("Dkappa:", dsctFcts.Dkappa)

dsctFcts.initA0Discrete()
print("a0:", dsctFcts.a0Discrete)

dsctFcts.initSigmaDiscrete()
print("Sigma:", dsctFcts.SigmaDiscrete)

dsctFcts.initSigma0Discrete()
print("sigma0:", dsctFcts.sigma0Discrete)

dsctFcts.initDSigma0Discrete()
print("Dsigma0:", dsctFcts.Dsigma0Discrete)


dsctFcts.initDKappa()
print("Dkappa:", dsctFcts.Dkappa)

# todo: number=10 -> complexWarning, unclear why.


wM = LinAl.WMatrix(dsctFcts)
print("fully initialized:", wM.dFfullyInitialized)

wM.initW0()
"""

def f(x):
	return x**2

if __name__ == '__main__':
	multiprocessing.freeze_support()
	p=multiprocessing.Pool(processes=4)
	result=p.apply_async(f,(2,))
	print(result.get())