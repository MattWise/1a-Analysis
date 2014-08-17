from __future__ import division

import math as m
from Utils import RadialCells

""" Representation of a set of initial conditions.
"""
class InitialParameters(object):
    def __init__(self, number=500, rStar=10 ** 11, rDisk=10 ** 15, p=1, q=5 / 8, mStar=10 ** 33, mDisk=10 ** 33,
                 a0Star=10**6, m=1):        # todo: getting a adequate value for a0Star
        self.number = number
        self.rStar = rStar
        self.rDisk = rDisk
        self.p = p
        self.q = q
        self.mStar = mStar
        self.mDisk = mDisk
        self.a0Star = a0Star
        self.G = 6.67384 * 10 ** (-11)
        self.m = m

    def __str__(self):
        return "number: " + str(self.number) \
               + "\n(rStar,rDisk): " + str((self.rStar, self.rDisk)) \
               + "\n(mStar,mDisk): " + str((self.mStar, self.mDisk)) \
               + "\n(p,q): " + str((self.p, self.q)) \
               + "\na0Star: " + str(self.a0Star) \
               + "\nm: " + str(self.m)


""" Class for calculating all the derived constants out there. Based on the initial constants.
"""
class DerivedConstants(object):
    def __init__(self, initParams):
        self.initialParameters = initParams
        self.n = 4 - (2 / initParams.q)
        self.sigmaStar = self.__sigStar(initParams)
        self.omegaStar = self.__omStar(initParams)
        self.qStar = self.__qStar(initParams)
        self.radialCells = self.__radCells(initParams)

    def __str__(self):
        return "n: " + str(self.n) + "\n" \
               + "sigmaStar: " + str(self.sigmaStar) + "\n" \
               + "omegaStar: " + str(self.omegaStar) + "\n" \
               + "qStar: " + str(self.qStar)

    def __sigStar(self, iP):
        return (2 - iP.p) * iP.mStar / \
               (2 * m.pi * iP.rStar ** 2 * ((iP.rDisk / iP.rStar) ** (2 - iP.p) - 1))

    def __omStar(self, iP):
        return m.sqrt(iP.G * iP.mStar / (iP.rStar ** 2))

    def __qStar(self, iP):
        return (self.omegaStar * iP.a0Star) / (m.pi * iP.G * self.sigmaStar)

    def __radCells(self, iP):
        return RadialCells(iP.number, iP.rStar, iP.rDisk)