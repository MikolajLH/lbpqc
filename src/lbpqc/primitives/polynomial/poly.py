from lbpqc.type_aliases import *

r'''
Polynomial is represented by a vector of it's coefficients, ordered by increasing powers.
E.g. X^3 + 7X - 2 = X^3 + 0X^2 + 7X + (-2) will be represented as:
[-2, 7, 0, 1]

Z[X] ring
'''




@enforce_type_check
def is_zero_poly(p: VectorInt) -> bool:
    r'''Checks if given polynomial is zero polynomial.

    Args:
        p: Polynomial's coefficients.
    
    Returns:
        True if $p \equiv 0$, False otherwise
    '''
    if len(p) == 0: raise ValueError("Empty numpy array is not a proper polynomial")
    
    return np.count_nonzero(p) == 0


@enforce_type_check
def deg(p: VectorInt,*, error_if_zero_poly: bool = False) -> int:
    r'''
    Returns degree of a given polynomial calculated as an index of the last nonzero ceofficient.
    Returns -1 as a degree of zero polynomial if `error_if_zero_poly` is set to `False`.

    Args:
        p: Polynomial's coefficients.
        error_if_zero_poly: Parameter deciding how to treat degree of zero polynomial.
    
    Returns:
        Degree of the given polynomial.
    
    Raises:
        ValueError: If given empty numpy array.
        ValueError: If given zero polynomial **and** error_if_zero_poly is set to True.
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
    so that resulting coefficient's arrray has length of $\deg(p) + 1$.

    If p is zero polynomial, then returns `np.array([0], dtype=int)`.

    Args:
        p: Polynomial's coefficients.
    
    Returns:
        Coefficient's array of length at most $\deg(p)$.
    
    Examples:
        >>> p = np.array([1,0,2,3,0,0,0])
        >>> trim(p)
        array([1,0,2,3])
    '''
    if is_zero_poly(p):
        return np.zeros(1, dtype=int)
    
    return p[:deg(p) + 1].copy()


@enforce_type_check
def pad(p: VectorInt, max_deg: int) -> VectorInt:
    r'''
    Pad's polynomial's coefficient's array with zero entries for powers higher than polynomial's degree,
    so that length of resulting array is equal to max_deg + 1.
    
    Args:
        p: Polynomial's coefficients.
        max_deg: Degree that $p$ is to be expanded to.
    
    Returns:
        Coefficient's array with length equal to `max_deg` + 1, filled with zeros at indices greater than $\deg(p)$.
    
    Examples:
        $p = X^3 + 7X - 2$
        >>> p = np.array([-2, 7, 0, 1])
        >>> pad(p, 5)
        array([-2, 7, 0, 1, 0, 0])
    '''
    if is_zero_poly(p):
        return zero_poly(max_deg)
    
    d = deg(p)
    if max_deg < d: raise ValueError("max_deg has to be greater or equal to the degree of a given polynomial p")
    
    return np.pad(trim(p), (0, max_deg - d))


@enforce_type_check
def monomial(coeff: int, degree: int) -> VectorInt:
    r'''
    For given degree $d$ and coefficient $c$, constructs a monomial
    $$
        cX^{d - 1}
    $$

    Args:
        coeff: Monomial's coefficient.
        degree: Monomial's degree.
    
    Returns:
        Coefficients' array with only nonzero entry `coeff` at `degree` index.
    
    Examples:
        $7X^5$
        >>> monomial(7, 5)
        array([0,0,0,0,0,5])
    '''
    p = np.zeros(degree + 1, dtype=int)
    p[degree] = coeff
    return p


@enforce_type_check
def zero_poly(max_deg: int = 0) -> VectorInt:
    r'''Explicitly constructs zero polynomial, i.e. a coefficient's array of length `max_deg` + 1 filled with zeros.

    Args:
        max_deg: .
    
    Returns:
        Coefficients' array of length `max_deg` + 1 filled with zeros.
    '''
    
    return np.zeros(max_deg + 1, dtype=int)



@enforce_type_check
def add(p: VectorInt, q: VectorInt) -> VectorInt:
    r'''Adds two polynomials.

    Args:
        p: polynomial's $p$ coefficients.
        q: polynomial's $q$ coefficients.
    
    Returns:
        Coefficients array of polynomial $p + q$.
    '''

    max_deg = max(deg(p), deg(q), 0)
    return trim(pad(p, max_deg) + pad(q, max_deg))


@enforce_type_check
def sub(p: VectorInt, q: VectorInt) -> VectorInt:
    r'''Subtract polynomial $q$ from polynomial $p$.

    Args:
        p: polynomial's $p$ coefficients.
        q: polynomial's $q$ coefficients.
    
    Returns:
        Coefficients array of polynomial $p - q$.
    '''

    max_deg = max(deg(p), deg(q), 0)
    return trim(pad(p, max_deg) - pad(q, max_deg))


@enforce_type_check
def mul(p: Vector, q: Vector) -> Vector:
    r'''Multiplies polynomials $p$ and $q$.
    
    Args:
        p: polynomial's $p$ coefficients.
        q: polynomial's $q$ coefficients.
    
    Returns:
        Coefficients array of polynomial $p \cdot q$.
    '''

    return np.polymul(p[::-1], q[::-1])[::-1]