""" This module covers the testing, debug and correctness implementations we did. Thus, it includes mathematical
(e.g. numerical) and physical testing functions encapsulated in classes needed for that purpose.
"""

import numpy as np

class Correctness(object):
	def __init__(self, eigenvalueSolver):
		self.eigenvalueSolver = eigenvalueSolver

	""" Returns a list of the (eigenvalue, eigvector) tuples that have regular
	eigenvalues.
	"""
	def regularIndices(self):
		if not(self.eigenvalueSolver.eigenvalues == None or self.eigenvalueSolver.eigenvectors == None):
			indicesOfRegularPairs = np.array([])
			for index in np.arange(len(self.eigenvalueSolver.eigenvalues)):
				iReal = self.eigenvalueSolver.eigenvalues[index].real
				iImag = self.eigenvalueSolver.eigenvalues[index].imag
				if not (iReal in [np.inf, -np.inf] \
						or iImag in [np.inf, -np.inf] \
						or iReal in [np.nan, -np.nan] \
						or iImag in [np.nan, -np.nan]):
					indicesOfRegularPairs = np.append(indicesOfRegularPairs, index)
			return indicesOfRegularPairs
		else:
			raise BaseException("eigenvalueSolver not initialized")


	""" Takes the given equation B14 and for each (eigenvalue, eigenvector)
	tuple whereas the eigenvalue is regular (e.g. no +/- inf and no +/- nan),
	the difference of both sides will be calculated. This should be close to
	zero if things work out.
	"""
	def testB14(self, verbose):
		if not(self.eigenvalueSolver.eigenvalues == None or self.eigenvalueSolver.eigenvectors == None):
			# todo: use np.array with correct working append
			results = []

			for index in self.regularIndices():
				matrix = self.eigenvalueSolver.wM.B14A \
				         - self.eigenvalueSolver.eigenvalues[index]*self.eigenvalueSolver.wM.B14C
				result = np.dot(matrix, self.eigenvalueSolver.eigenvectors[index])
				results.append(result)

			if verbose:
				abs = map(np.linalg.norm, results)
				print "testB14(): minimum abs =", min(abs)
				print "testB14(): maximum abs =", max(abs)

			return results

		else:
			raise BaseException("eigenvalueSolver not initialized")

	""" Takes the calculated eigenvalues and -vectors in order to test if they satisfy equation
	B12.
	"""
	def testB12(self, verbose):
		if not(self.eigenvalueSolver.eigenvalues == None or self.eigenvalueSolver.eigenvectors == None):

			number = self.eigenvalueSolver.wM.number
			W0 = self.eigenvalueSolver.wM.W0
			W1 = self.eigenvalueSolver.wM.W1
			W2 = self.eigenvalueSolver.wM.W2
			W3 = self.eigenvalueSolver.wM.W3
			W4 = self.eigenvalueSolver.wM.W4
			W5 = self.eigenvalueSolver.wM.W5

			results = []

			for index in self.regularIndices():
				omega = self.eigenvalueSolver.eigenvalues[index]
				matrix = W0 + omega * W1 \
				         + omega**2 * W2 \
				         + omega**3 + W3 \
				         + omega**4 + W4 \
				         + omega**5 * W5
				result = np.dot(matrix, self.eigenvalueSolver.eigenvectors[index][0:number])
				results.append(result)

			if verbose:
				abs = map(np.linalg.norm, results)
				print "testB12(): minimum abs =", min(abs)
				print "testB12(): maximum abs =", max(abs)

			return results

		else:
			raise BaseException("eigenvalueSolver not initialized")