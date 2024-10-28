import numpy as np

r'''
Hint for elements of $\\mathbb{Z}_{p}$ for some modulus $q$.
'''
type ModInt = int


r'''
Hints for elements of $\\mathbb{R}^{n}$ vector space.
Components types are not explicitly defined.
'''
type Vector = np.ndarray
type Matrix = np.ndarray
type SquareMatrix = np.ndarray


r'''
Hints for elements of $\\mathbb{R}^{n}$ vector space.
Components types are are explicitly floats.
'''
type VectorFloat = np.ndarray[float]
type MatrixFloat = np.ndarray[float]
type SquareMatrixFloat = np.ndarray[float]


r'''
Hints for elements of $\\mathbb{Z}^{n}$ vector space.
'''
type VectorInt = np.ndarray[int]
type MatrixInt = np.ndarray[int]
type SquareMatrixInt = np.ndarray[int]


r'''
Hints for elements of $\\mathbb{Z}_{p}^{n}$ vector space for some modulus $q$.
'''
type VectorMod = np.ndarray[int]
type MatrixMod = np.ndarray[int]
type SquareMatrixMod = np.ndarray[int]
