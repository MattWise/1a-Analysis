from __future__ import division

import Functions, LinAl, Params

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

dsctFcts.initA0Discrete()
print "a0:", dsctFcts.a0Discrete

dsctFcts.initSigmaDiscrete()
print "Sigma:", dsctFcts.SigmaDiscrete

dsctFcts.initSigma0Discrete()
print "sigma0:", dsctFcts.sigma0Discrete

dsctFcts.initDSigma0Discrete()
print "Dsigma0:", dsctFcts.Dsigma0Discrete

dsctFcts.initDKappa()
print "Dkappa:", dsctFcts.Dkappa

# todo: number=10 -> complexWarning, unclear why.

wM = LinAl.WMatrix(dsctFcts)
print "fully initialized:", wM.dFfullyInitialized

wM.initW0()
