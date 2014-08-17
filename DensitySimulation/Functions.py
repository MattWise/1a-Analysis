from __future__ import division

import math as m
import numpy as np
import Utils as U

""" This class contains the functions we have in analytical form. Contains non-static methods only because
of dependency on the initial parameters.
"""
class AnalyticalFunctions(object):
    def __init__(self, iP, dC):
        self.iP = iP
        self.dC = dC

    # unpurtubed disk surface density
    def sigma0(self,r):
        return self.dC.sigmaStar*(self.iP.rStar/float(r))**self.iP.p

    # sound speed
    def a0(self,r):
        return self.iP.a0Star*float(r)**(-self.iP.q/2)

    # big sigma
    def Sigma(self,r):
        return (2*m.pi*self.iP.G*self.sigma0(float(r)))/(self.a0(float(r))**2)


""" This class represents everything we have in a discrete manner. Besides that, a - also public available -
discretizer is available based on the given parameters.
"""
class DiscreteFunctions(object):
    def __init__(self, analyticFcts):
        self.iP = analyticFcts.iP
        self.dC = analyticFcts.dC
        self.analyticFunctions = analyticFcts
        self.a0Discrete = np.array([])
        self.sigmaDiscrete = np.array([])
        self.sigma0Discrete = np.array([])
        self.Omega = np.array([])
        self.kappa = np.array([])
        self.linAlg = U.LinearAlgebra(self.iP, self. dC) # todo: possible reason to do it like that is: default
        # todo: maths can be used by user without taking a look at it but he can also change it!
        self.linAlg.initLDMOne()
        self.linAlg.initLDMTwo()

    # =============== discretizing function ===============
    def discretizer(self, function):
        res = np.array([])
        for cellR in self.dC.radialCells.rValues:
            res = np.append(res, function(cellR))
        return res

    # ================ discrete functions =================
    # todo: dummy version
    def calcOmega(self):
        return np.ones(self.iP.number)

    def calcKappa(self):
        prefactor = (np.asarray(dC.radialCells.rValues))**(-4.0)
        toDerive = np.asmatrix(((np.asarray(dC.radialCells.rValues)**2)*self.Omega)**2).T
        return (prefactor * (self.linAlg.logDerivationMatrixOne * toDerive))
        # todo: does not work

    # ============== init discretized fields ==============
    def initA0Discrete(self):
        if len(self.a0Discrete) == 0:
            self.a0Discrete = self.discretizer(self.analyticFunctions.a0)
        else:
            raise BaseException("already initialized")

    def initSigmaDiscrete(self):
        if len(self.sigmaDiscrete) == 0:
            self.sigmaDiscrete = self.discretizer(self.analyticFunctions.Sigma)
        else:
            raise BaseException("already initialized")

    def initSigma0Discrete(self):
        if len(self.sigma0Discrete) == 0:
            self.sigma0Discrete = self.discretizer(self.analyticFunctions.sigma0)
        else:
            raise BaseException("already initialized")

    def initOmega(self):
        if len(self.Omega) == 0:
            self.Omega = self.calcOmega()
        else:
            raise BaseException("already initialized")

    def initKappa(self):
        if len(self.kappa) == 0:
            self.kappa = self.calcKappa()
        else:
            raise BaseException("already initialized")





# ======================== TESTING AREA ========================

import Params as P
iP = P.InitialParameters(number = 10)
dC = P.DerivedConstants(iP)
anF = AnalyticalFunctions(iP, dC)
diF = DiscreteFunctions(anF)

diF.initOmega()
diF.initA0Discrete()
diF.initKappa()

print diF.kappa

