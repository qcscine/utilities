__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as scine
import numpy as np
import os

class SigmaVectorEvaluatorPython(scine.SigmaVectorEvaluator):
    def __init__(self, matrix):
        scine.SigmaVectorEvaluator.__init__(self)
        self.matrix = matrix
    def evaluate(self, guess_vectors):
        return np.dot(self.matrix, guess_vectors)
    def collapsed(self, newSubspaceDimension):
        return
    def swap(self, i, j):
        return

def create_matrix():
    # create a selfadjoint matrix
    matrix = np.random.rand(100,100)
    matrix = 0.5*(matrix + np.transpose(matrix))
    matrix[np.diag_indices_from(matrix)] += 1
    return matrix

def initialize_diagonalizer(matrix):
    # Create sigma vector evaluator and preconditioner
    sve = scine.IndirectSigmaVectorEvaluator(matrix)
    prec = scine.IndirectPreconditionerEvaluator(matrix[np.diag_indices_from(matrix)])

    # Create and fill Non Orthogonal Davidson
    diag = scine.NonOrthogonalDavidson(5,100)
    diag.sigma_vector_evaluator = sve
    diag.set_preconditioner(prec)

    return diag

def test_SigmaVectorEvaluator():
    ref = create_matrix()
    sve = scine.IndirectSigmaVectorEvaluator(ref)
    result = sve.evaluate(2.0 * np.identity(100))
    assert np.all(2.0 * ref[:,:] == result[:,:])

def test_Preconditioner():
    '''
    Test that if you try to precondition a vector of ones, you just get
    -1.0 / (difference btw the diagonal and the current eigenvalue)
    '''
    ref = create_matrix()
    diag = ref[np.diag_indices_from(ref)]
    ones_vector = np.ones(100)
    arbitrary_eigenvalue = 3.5
    prec = scine.IndirectPreconditionerEvaluator(diag)
    result = prec.evaluate(ones_vector, arbitrary_eigenvalue)
    assert np.all(result[:] == -1.0 / (diag - arbitrary_eigenvalue))

def test_InitializeDiagonalizer():
    diag = initialize_diagonalizer(create_matrix())

def test_DiagonalizeWithNonOrthogonalDavidson():
    ref = create_matrix()
    diag = initialize_diagonalizer(ref)
    result = diag.solve(scine.core.Log.silent())
    # Get reference numbers
    w, v = np.linalg.eig(ref)
    assert np.all(result.eigenvalues[:] - sorted(w)[:5] <= 1.0e-5)

def test_DiagonalizeWithOrthogonalDavidson():
    ref = create_matrix()
    # Create sigma vector evaluator and preconditioner
    sve = scine.IndirectSigmaVectorEvaluator(ref)
    prec = scine.IndirectPreconditionerEvaluator(ref[np.diag_indices_from(ref)])

    # Create and fill Non Orthogonal Davidson
    diag = scine.OrthogonalDavidson(5,100)
    diag.sigma_vector_evaluator = sve
    diag.set_preconditioner(prec)
    result = diag.solve(scine.core.Log.silent())
    # Get reference numbers
    w, v = np.linalg.eig(ref)
    assert np.all(result.eigenvalues[:] - sorted(w)[:5] <= 1.0e-5)

def test_DiagonalizeWithPythonSigmaVectorEvaluator():
    ref = create_matrix()
    diag = initialize_diagonalizer(ref)
    # Set python specific sigma vector evaluator
    # Note: first initialize, then assign to prevent auto casting.
    #       If I write diag.sigma_vector_evaluator = SigmaVectorEvaluatorPython(ref)
    #       then it it tried to look for the method SigmaVectorEvaluator::evaluate()
    #       instead of SigmaVectorEvaluatorPython::evaluate()
    sve = SigmaVectorEvaluatorPython(ref)
    diag.sigma_vector_evaluator = sve
    result = diag.solve(scine.core.Log.silent())
    # Get reference numbers
    w, v = np.linalg.eig(ref)
    assert np.all(result.eigenvalues[:] - sorted(w)[:5] <= 1.0e-5)
