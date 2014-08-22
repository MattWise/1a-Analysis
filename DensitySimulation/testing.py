from __future__ import division,print_function
import Functions, LinAl, Params
import numpy as np

def logDerivationMatrix(dimension, order, innerBoundary=False, outerBoundary=False):

	#creates a differentiation matrix for numerical differentiation. "Order=0" means identity matrix, "order=1" means first derivative, "order=2" means second.


	def diagonals(array,a,b,value):
		#goes down diagonal form specified coordinates, setting each to value.
		if a>b:
			times = dimension-a-1
		else:
			times=dimension-b
		for row in range(times):
			array[a,b]=value
			a+=1
			b+=1
		return array


	def setBoundary(array,dimension,innerBoundary,outerBoundary):
		#sets the matrix boundary conditions. Defaults are from ARS paper.

		outerBoundary.reverse()

		for index,value in enumerate(innerBoundary):
			array[0,index]=value
		for index,value in enumerate(outerBoundary):
			array[dimension-1,dimension-1-index]=value
		return array

	array = np.zeros([dimension,dimension])
	prefactor = 0

	if order==0:
		array=diagonals(array,0,0,1)
		prefactor=np.complex128(1)

	if order == 1:
		array = diagonals(array,1,0,-1)
		array = diagonals(array,1,2,1)

		#defualt boundary conditions
		if not innerBoundary:
			inner_boundary=[-3,4,-1]
		if not outerBoundary:
			outer_boundary=[1,-4,3]

		array=setBoundary(array,dimension,inner_boundary,outer_boundary)
		prefactor = np.complex128(1/(2*np.log10(self.derivedConstants.radialCells.f)))

	if order == 2:
		array = diagonals(array,1,0,1)
		array = diagonals(array,1,1,-2)
		array = diagonals(array,1,2,1)

		#defualt boundary conditions
		if not innerBoundary:
			inner_boundary = [2,-5,4,-1]
		if not outerBoundary:
			outer_boundary=[-14,-5,2]

		array=setBoundary(array,dimension,inner_boundary,outer_boundary)
		prefactor = np.complex128(1/(np.log10(self.derivedConstants.radialCells.f))**2)


	return prefactor * array

print(logDerivationMatrix(100,0))

"""
iP = Params.InitialParameters(number=100)
print str(iP)
dC = Params.DerivedConstants(iP)
print str(dC)

anaFcts = Functions.AnalyticalFunctions(iP, dC)




dsctFcts = Functions.DiscreteFunctions(anaFcts)

dsctFcts.initOmega()
print "Omega:", dsctFcts.Omega

dsctFcts.initDOmega()
print "DOmega:", dsctFcts.DOmega

dsctFcts.initDDOmega()
print "DDOmega:", dsctFcts.DDOmega

dsctFcts.initKappa()
print "kappa:", dsctFcts.kappa

dsctFcts.initDkappa()
print "Dkappa:", dsctFcts.Dkappa

dsctFcts.initA0Discrete()
print "a0:", dsctFcts.a0Discrete

dsctFcts.initSigmaDiscrete()
print "Sigma:", dsctFcts.SigmaDiscrete

dsctFcts.initSigma0Discrete()
print "sigma0:", dsctFcts.sigma0Discrete

dsctFcts.initDSigma0Discrete()
print "Dsigma0:", dsctFcts.Dsigma0Discrete

# ===== works up to here =====










# ===== what does not work follows: ======



"""