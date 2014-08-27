from __future__ import division,print_function

import Functions, LinAl, Params
import numpy as np
import multiprocessing
import time as t

import matplotlib.pyplot as plt

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

eigSol = LinAl.EigenvalueSolver(wM)



eigSol.initEigen()

print("# of eigenvalues", len(eigSol.eigenvalues))

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

# plt.scatter(real[4:],img[4:])
plt.scatter(real, img)
plt.show()