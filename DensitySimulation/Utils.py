from __future__ import division

import numpy as np
import time as t
import math as m


""" This class represents a set of radial grid points having given
the number of cells and radius of star & radius of disk. Uses numpy arrays.
"""
class RadialCells(object):
    def __init__(self, number, rStar, rDisk):
        # simple constructor without setters, getters and stuff. thus, the user is
        # is responsible to use the class in a proper way.
        self.number = float(number)
        self.rStar = float(rStar)
        self.rDisk = float(rDisk)
        self.f = self.__calcF()
        self.rValues = None
        self.__calcRvalues()

    # ===============
    # inner methods:
    # ===============

    # definition see page 973, B1b
    def __calcF(self):
        return (self.rDisk / self.rStar) ** (1 / (self.number - 1))

    # definition see page 972, B1a
    def __calcRvalues(self):
        self.rValues = np.array([])
        for i in range(1, int(self.number)):
            value = self.f ** (i - 1) * self.rStar
            self.rValues = np.append(self.rValues, value)
        # here, we assert that r_N = rDisk for completness but
        # due to numerical errors, we have to fix that the following
        # way
        self.rValues = np.append(self.rValues, self.rDisk)

    # simple toString method for proper debug and information
    def __str__(self):
        return "number: " + str(self.number) + "\n" \
               + "rStar: " + str(self.rStar) + "\n" \
               + "rDisk: " + str(self.rDisk)


""" This class holds all the matrices needed for numerical derivatives and integration. Since matrix creation is
quite intense in time and power, we use a non-static class here in order to store the calculated matrices.
We use Numpy.Matrix for performance. The class needs to be initialized in order to inform the user that some
code might take some time - no exception handling here only for multiple-inits.
"""
class LinearAlgebra:
    def __init__(self, iP, dC):
        # basic parameters
        self.initialParameters = iP
        self.derivedConstants = dC
        # first, not initialized to save time. besides that,
        # these are the checking values if already initialized.
        # but there are no control mechanisms.
        self.logDerivationMatrixOne = None
        self.logDerivationMatrixTwo = None

    # ======== init methods ========

    """ Initialization of logarithmic matrix of order 1.
    """
    def initLDMOne(self):
        if self.logDerivationMatrixOne == None:
            self.logDerivationMatrixOne\
                = self.__logDerivationMatrix(self.initialParameters.number, 1) # use method with correct input.
        else:
            raise BaseException("logDerivationMatrix of order=1 already initialized")

    """ Initialization of logarithmic matrix of order 2.
    """
    def initLDMTwo(self):
        if self.logDerivationMatrixTwo == None:
            self.logDerivationMatrixTwo\
                = self.__logDerivationMatrix(self.initialParameters.number, 2) # use method with correct input.
        else:
            raise BaseException("logDerivationMatrix of order=2 already initialized")

    """ For a given dimension size N, returns the logarithmic derivation matrices
    of order one and two: the discretized version of d/dlog(r).
    used to get d/dr and (d/dr)^2, derivation see page 973, B3a and B3b.
    """
    def __logDerivationMatrix(self, dimension, order):
        # TODO: matts method as replacement because faster.
        # first, create a numpy.array with all the numbers and
        # then, just convert it to numpy.matrix. therefore, create
        # a two dimensional list.
        if not (order == 1 or order == 2):
            raise BaseException("order must be 1 or 2")
        if dimension < 4:
            raise BaseException("dimension must be greater 3")
        res = []
        prefactor = {1: 1/(2*m.log10(self.derivedConstants.radialCells.f)),
                     2: (1/(2*m.log10(self.derivedConstants.radialCells.f)))**2}.get(order)
        for i in range(dimension):
            # helping switch case functions
            def init1(n):
                return {1: [-3, 4, -1], 2: [2, -5, 4, -1]}.get(n)
            def init2rev(n):
                return {1: [3, -4, 1], 2: [2, -5, -14]}.get(n)
            def nrs(n):
                return {1: [-1, 0, 1], 2: [1, -2, 1]}.get(n)
            # init
            initial1 = init1(order)
            initial2reverse = init2rev(order)
            numbers = nrs(order)
            row = []
            # filling
            if i == 0:
                row.extend(initial1)
                for j in range(dimension - len(row)):
                    row.append(0)
            elif i == dimension - 1:
                row.extend(initial2reverse)
                for j in range(dimension - 3):
                    row.append(0)
                row.reverse()
            else:
                for j1 in range(0, i - 1):
                    row.append(0)
                row.extend(numbers)
            for j2 in range(len(row), dimension):
                row.append(0)
            res.append(row)
        return prefactor * np.matrix(res)
