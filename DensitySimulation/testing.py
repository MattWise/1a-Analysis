from __future__ import division,print_function

import Functions, LinAl, Params, Correctness
import numpy as np
import multiprocessing
import time as t
import sys

import matplotlib.pyplot as plt

# constants
seperator = "========================================================\n\n\n========================================================"

# numpy settings
np.set_printoptions(threshold=sys.maxint, edgeitems=sys.maxint, linewidth=sys.maxint)
np.seterr(all="ignore") # this setting reveals the infinity problem

def timing(function):
	timeStart = t.clock()
	function()
	timeEnd = t.clock()
	print("took me:", (timeEnd-timeStart), "s")

def determinants(wMatrix):
	tmp = [wMatrix.W0, wMatrix.W1, wMatrix.W2, wMatrix.W3, wMatrix.W4, wMatrix.W5]
	for i in range(len(tmp)):
		print("det(W" + str(i) + ")=", str(np.linalg.det(tmp[i])))
	print("det(B14A)=", str(np.linalg.det(wMatrix.B14A)))
	print("det(B14C)=", str(np.linalg.det(wMatrix.B14C)))

def normalVSbadEigenvalues(eigenvalues):
	normal=0
	bad=0
	for i in eigenvalues:
		iReal = i.real
		iImag = i.imag
		if iReal in [np.inf, -np.inf] \
				or iImag in [np.inf, -np.inf] \
				or iReal in [np.nan, -np.nan] \
				or iImag in [np.nan, -np.nan]:
			bad += 1
		else:
			normal += 1
	print(normal, bad)

iP = Params.InitialParameters(number=100)
dC = Params.DerivedConstants(iP)


anaFcts = Functions.AnalyticalFunctions(iP, dC)
dsctFcts = Functions.DiscreteFunctions(anaFcts)

dsctFcts.init()

wM = LinAl.WMatrix(dsctFcts)

def init():
	wM.initW0()
	wM.initW1()
	wM.initW2()
	wM.initW3()
	wM.initW4()
	wM.initW5()
	wM.initB14A()
	wM.initB14C()

init()
print("wM init() finished")

print(seperator)

#determinants(wM)


eigSol = LinAl.EigenvalueSolver(wM)
eigSol.initEigen()

correctness = Correctness.Correctness(eigSol)

correctness.testB12(True)
correctness.testB14(True)




"""
real = []
img = []

for z in eigSol.eigenvalues:
	rp = z.real
	im = z.imag
	if not (rp == np.inf or rp == - np.inf) \
			and not (im == np.inf or im == - np.inf):
		real.append(rp)
		img.append(im)


print("real",real)
print("realLen=" + str(len(real)))
print("img",img)
print("imagLen=" + str(len(img)))

print(min(real), max(real))
print(min(img), max(img))

"""
"""
# plt.scatter(real[4:],img[4:])
plt.scatter(real, img)
plt.grid(True)
plt.xlabel("realpart")
plt.ylabel("imagpart")
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.show()
"""