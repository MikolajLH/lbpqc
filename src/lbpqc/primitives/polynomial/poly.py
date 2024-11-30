from lbpqc.type_aliases import *
from typing import Tuple

r'''
Polynomial is represented by a vector of it's coefficients, ordered by increasing powers.
E.g. X^3 + 7X - 2 = X^3 + 0X^2 + 7X + (-2) will be represented as:
[-2, 7, 0, 1]

Z[X] ring
'''





def is_zero_poly(p: VectorInt) -> bool:
    r'''
    Checks if given polynomial is zero polynomial
    '''
    if not is_VectorInt(p): raise TypeError()
    if len(p) == 0: raise ValueError("Empty numpy array is not a proper polynomial")
    
    return np.count_nonzero(p) == 0


def deg(p: VectorInt,*, error_if_zero_poly: bool = False) -> int:
    r'''
    Returns degree of a given polynomial calculated as an index of the last nonzero ceofficient.
    If error_if_zero_poly is set to True, raises ValueError, 
    else returns -1 as a degree of zero polynomial.
    Raises ValueError if given an empty numpy array
    '''
    if not is_VectorInt(p): raise TypeError()
    
    if len(p) == 0: raise ValueError("degree undefined for an empty numpy array")
    
    if len(nonzeros := np.nonzero(p)[0]) == 0:
        if error_if_zero_poly: raise ValueError("Degree of zero polynomial is undefined")
        return -1
    else:
        return nonzeros[-1]



def trim(p: VectorInt) -> VectorInt:
    r'''
    Trims zero coefficients of powers higher than polynomial's degree,
    so that len(res) == deg(p) + 1.
    If p is zero polynomial returns np.array([0], dtype=int).
    '''
    if not is_VectorInt(p): raise TypeError()

    if is_zero_poly(p):
        return np.zeros(1, dtype=int)
    
    return p[:deg(p) + 1].copy()


def pad(p: VectorInt, max_deg: int) -> VectorInt:
    r'''
    adds zero coefficients of powers higher than polynomial's degree,
    so that len(output) == max_deg + 1
    E.g.
    p = X^3 + 7X - 2 -> p := [-2, 7, 0, 1]
    pad(p, 5) -> [-2, 7, 0 ,1, 0, 0] -> 0X^5 + 0X^4 + X^3 + 0X^2 + 7X - 2
    '''
    if not is_VectorInt(p): raise TypeError()
    if is_zero_poly(p):
        return zero_poly(max_deg)
    
    d = deg(p)
    if max_deg < d: raise ValueError("max_deg has to be greater or equal to the degree of a given polynomial p")
    
    return np.pad(trim(p), (0, max_deg - d))


def monomial(coeff: int, degree: int) -> VectorInt:
    r'''
    Constructs a monomial of a given degree with a given coefficient
    '''
    p = np.zeros(degree + 1, dtype=int)
    p[degree] = coeff
    return p


def zero_poly(max_deg: int = 0) -> VectorInt:
    return np.zeros(max_deg + 1, dtype=int)


def add(p: VectorInt, q: VectorInt) -> VectorInt:
    if not is_VectorInt(p) or not is_VectorInt(q): raise TypeError()
    
    max_deg = max(deg(p), deg(q), 0)
    return trim(pad(p, max_deg) + pad(q, max_deg))


def sub(p: VectorInt, q: VectorInt) -> VectorInt:
    if not is_VectorInt(p) or not is_VectorInt(q): raise TypeError()
    
    max_deg = max(deg(p), deg(q), 0)
    return trim(pad(p, max_deg) - pad(q, max_deg))


def mul(p: Vector, q: Vector) -> Vector:
    if not is_VectorInt(p) or not is_VectorInt(q): raise TypeError()

    return np.polymul(p[::-1], q[::-1])[::-1]