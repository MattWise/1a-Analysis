from __future__ import division,print_function

import Functions, LinAl, Params
import numpy as np
import time

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

def calcB14C():
	zeros=np.zeros((500,500),dtype=np.complex128)
	I=np.identity(500,dtype=np.complex128)
	C2=np.concatenate((zeros,zeros,I,zeros,zeros),0)
	C3=np.concatenate((zeros,zeros,zeros,I,zeros),0)
	C4=np.concatenate((zeros,zeros,zeros,zeros,I),0)
	C5=np.concatenate((I,zeros,zeros,zeros,zeros),0)
	columns=(C2,C3,C4,C5)
	print()
	return np.concatenate(columns,1)

A=time.clock()
print(calcB14C())
B=time.clock()
print(B-A)