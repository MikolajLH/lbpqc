import numpy as np
from typing import Any



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



r'''
Predicates for type checking
'''
def is_nparray(obj: Any) -> bool:
    return isinstance(obj, np.ndarray)


def is_Vector(obj: Any) -> bool:
    return is_nparray(obj) and len(obj.shape) == 1


def is_Matrix(obj: Any) -> bool:
    return is_nparray(obj) and len(obj.shape) == 2


def is_SquareMatrix(obj: Any) -> bool:
    return is_Matrix(obj) and obj.shape[0] == obj.shape[1]


def is_VectorInt(obj: Any) -> bool:
    return is_Vector(obj) and obj.dtype == int


def is_MatrixInt(obj: Any) -> bool:
    return is_Matrix(obj) and obj.dtype == int


def is_SquareMatrixInt(obj: Any) -> bool:
    return is_SquareMatrix(obj) and obj.dtype == int



def is_VectorFloat(obj: Any) -> bool:
    return is_Vector(obj) and obj.dtype == float


def is_MatrixFloat(obj: Any) -> bool:
    return is_Matrix(obj) and obj.dtype == float


def is_SquareMatrixFloat(obj: Any) -> bool:
    return is_SquareMatrix(obj) and obj.dtype == float
