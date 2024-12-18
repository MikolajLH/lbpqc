from lbpqc.type_aliases import *

r'''
Polynomial is represented by a vector of it's coefficients, ordered by increasing powers.
E.g. X^3 + 7X - 2 = X^3 + 0X^2 + 7X + (-2) will be represented as:
[-2, 7, 0, 1]

Z[X] ring
'''




@enforce_type_check
def is_zero_poly(p: VectorInt) -> bool:
    r'''
    Checks if given polynomial is zero polynomial
    '''
    if len(p) == 0: raise ValueError("Empty numpy array is not a proper polynomial")
    
    return np.count_nonzero(p) == 0


@enforce_type_check
def deg(p: VectorInt,*, error_if_zero_poly: bool = False) -> int:
    r'''
    Returns degree of a given polynomial calculated as an index of the last nonzero ceofficient.
    If error_if_zero_poly is set to True, raises ValueError, 
    else returns -1 as a degree of zero polynomial.
    Raises ValueError if given an empty numpy array
    '''
    if len(p) == 0: raise ValueError("degree undefined for an empty numpy array")
    
    if len(nonzeros := np.nonzero(p)[0]) == 0:
        if error_if_zero_poly: raise ValueError("Degree of zero polynomial is undefined")
        return -1
    else:
        return nonzeros[-1]


@enforce_type_check
def trim(p: VectorInt) -> VectorInt:
    r'''
    Trims zero coefficients of powers higher than polynomial's degree,
    so that len(res) == deg(p) + 1.
    If p is zero polynomial returns np.array([0], dtype=int).
    '''
    if is_zero_poly(p):
        return np.zeros(1, dtype=int)
    
    return p[:deg(p) + 1].copy()


@enforce_type_check
def pad(p: VectorInt, max_deg: int) -> VectorInt:
    r'''
    adds zero coefficients of powers higher than polynomial's degree,
    so that len(output) == max_deg + 1
    E.g.
    p = X^3 + 7X - 2 -> p := [-2, 7, 0, 1]
    pad(p, 5) -> [-2, 7, 0 ,1, 0, 0] -> 0X^5 + 0X^4 + X^3 + 0X^2 + 7X - 2
    '''
    if is_zero_poly(p):
        return zero_poly(max_deg)
    
    d = deg(p)
    if max_deg < d: raise ValueError("max_deg has to be greater or equal to the degree of a given polynomial p")
    
    return np.pad(trim(p), (0, max_deg - d))


@enforce_type_check
def monomial(coeff: int, degree: int) -> VectorInt:
    r'''
    Constructs a monomial of a given degree with a given coefficient
    '''
    p = np.zeros(degree + 1, dtype=int)
    p[degree] = coeff
    return p


@enforce_type_check
def zero_poly(max_deg: int = 0) -> VectorInt:
    
    return np.zeros(max_deg + 1, dtype=int)



@enforce_type_check
def add(p: VectorInt, q: VectorInt) -> VectorInt:

    max_deg = max(deg(p), deg(q), 0)
    return trim(pad(p, max_deg) + pad(q, max_deg))


@enforce_type_check
def sub(p: VectorInt, q: VectorInt) -> VectorInt:

    max_deg = max(deg(p), deg(q), 0)
    return trim(pad(p, max_deg) - pad(q, max_deg))


@enforce_type_check
def mul(p: Vector, q: Vector) -> Vector:

    return np.polymul(p[::-1], q[::-1])[::-1]