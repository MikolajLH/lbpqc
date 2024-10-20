import numpy as np
from numpy.polynomial import polynomial as np_poly
from lbpqc.primitives.integer import integer_ring


def has_int_coeffs(coeffs: np.ndarray[int]) -> bool:
    r'''
    Predicate used in assertions. Checks if numpy array's contained type is int.
    '''
    return coeffs.dtype == int

def is_zero_poly(coeffs: np.ndarray[int], modulus: int) -> bool:
    r'''
    Checks if the polynomial with given coefficients is a zero polynomial in a ring.
    '''
    assert has_int_coeffs(coeffs)

    return not np.any(reduce_coeffs_modulo(coeffs, modulus))



def is_const_poly(coeffs: np.ndarray) -> bool:
    r'''
    Checks if the polynomial $ p(X) $ with given coefficients is a constant polynomial $ p(X) \equiv C \quad C \in \mathcal{Z}\\{0\} $.
    '''
    nz = np.nonzero(coeffs)[0]
    # there is only one nonzero coefficient and it's constant coefficient
    return len(nz) == 1 and nz[0] == 0



def ring_zero() -> np.ndarray[int]:
    r'''
    Returns zero polynomial in a ring.
    '''
    return np.array([0])


def poly_monomial(m_deg: int, coeff: int) -> np.ndarray[int]:
    r'''
    Constructs and returns monomial $m(X)$ of a given degree $d$ with a given coefficient $c$.

    $$
    m(X) = c \cdot X^{d}
    $$
    '''
    p = np.zeros(m_deg + 1, dtype=int)
    p[m_deg] = coeff
    return p



def poly_NTRU_modulus(N: int) -> np.ndarray[int]:
    r'''
    Construct and returns polynomial $p(X)$ used in NTRU cryptosystem.

    $$
    p(X) = X^{N} - 1
    $$
    '''
    p = np.zeros(N + 1, dtype=int)
    p[0] =-1
    p[N] = 1
    return p


def deg(f: np.ndarray[int], modulus: int) -> int:
    r'''
    Returns index of last nonzero element of a vector.

    Degree of 0 polynomial is undefined.
    '''
    assert has_int_coeffs(f)
    assert not is_zero_poly(f, modulus)

    return np.nonzero(f)[0][-1]

def trim_trailing_zeros(p: np.ndarray[int], modulus: int) -> np.ndarray[int]:
    r'''
    Omits coefficients of powers higher than degree of given polynomial.

    Returns a vector which dimension is equal to its degree.
    '''
    assert has_int_coeffs(p)

    if is_zero_poly(p, modulus):
        return ring_zero()
    
    return p[:deg(p, modulus) + 1]


def pad_with_trailing_zeros(p: np.ndarray[int], n: int) -> np.ndarray[int]:
    r'''
    Returns vector which dimension is equal to parameter $n$ regardless of the degree of polynomial.

    Adds trailing 0 coefficients for powers higher than the degree.
    '''
    assert has_int_coeffs(p)
    assert n >= len(p)

    return np.pad(p, (0, n - len(p)))


def reduce_coeffs_modulo(p: np.ndarray[int], modulus: int) -> np.ndarray[int]:
    r'''
    Reduce coefficients modulo **modulus** so $ 0 \leq a_i < modulus $.
    '''
    assert has_int_coeffs(p)

    return p % modulus


def poly_mod_add(f: np.ndarray[int], g: np.ndarray[int], modulus: int) -> np.ndarray[int]:
    r'''
    Addition of two polynomial in ring $\mathcal{Z}_{p}[X]$.
    '''
    assert has_int_coeffs(f)
    assert has_int_coeffs(g)

    return reduce_coeffs_modulo(np_poly.polyadd(f, g).astype(int), modulus)

def poly_mod_sub(f: np.ndarray[int], g: np.ndarray[int], modulus: int) -> np.ndarray[int]:
    r'''
    Substraction of two polynomial in ring $\mathcal{Z}_{p}[X]$.
    '''
    assert has_int_coeffs(f)
    assert has_int_coeffs(g)

    return reduce_coeffs_modulo(np_poly.polysub(f, g).astype(int), modulus)


def poly_mod_mul(f: np.ndarray[int], g: np.ndarray[int], modulus: int) -> np.ndarray[int]:
    r'''
    Multiplication of two polynomial in ring $\mathcal{Z}_{p}[X]$.
    '''
    assert has_int_coeffs(f)
    assert has_int_coeffs(g)

    return reduce_coeffs_modulo(np_poly.polymul(f, g).astype(int), modulus)

def poly_mod_scalar_mul(p: np.ndarray[int], s: int, modulus: int) -> np.ndarray[int]:
    r'''
    Multiplication by scalar ring $\mathcal{Z}_{p}[X]$.
    '''
    assert has_int_coeffs(p)

    return reduce_coeffs_modulo(s * p, modulus)


def poly_mod_euclidean_division(f: np.ndarray[int], g: np.ndarray[int], modulus) -> tuple[np.ndarray[int], np.ndarray[int]]:
    r'''
    Implementation of Euclidean division algorithm in ring $\mathcal{Z}_{p}[X]$.
    '''
    assert has_int_coeffs(f)
    assert has_int_coeffs(g)
    assert not is_zero_poly(g, modulus)

    q = ring_zero()
    r = f.copy()
    while (not is_zero_poly(r, modulus)) and deg(r, modulus) >= deg(g, modulus):
        n, m = deg(r, modulus), deg(g, modulus)
        gm_inv = integer_ring.modular_inverse(g[m], modulus)
        c = (r[n] * gm_inv) % modulus

        q = poly_mod_add(q, poly_monomial(n - m, c), modulus)
        r = poly_mod_sub(r, poly_mod_mul(g, poly_monomial(n - m, c), modulus), modulus)

    return q, r


def poly_mod_gcd(f: np.ndarray[int], g: np.ndarray[int], modulus: int) -> np.ndarray[int]:
    r'''
    Computes gcd of two polynomials in ring $\mathcal{Z}_{p}[X]$.
    '''
    assert has_int_coeffs(f)
    assert has_int_coeffs(g)

    if deg(g, modulus) > deg(f, modulus):
        f, g = g, f
    while not is_zero_poly(g, modulus):
        _, r =  poly_mod_euclidean_division(f, g, modulus)
        f = g
        g = r
    return f


def poly_mod_coprime(f: np.ndarray[int], g: np.ndarray[int], modulus: int) -> bool:
    r'''
    Check if two polynomials are coprime in ring $\mathcal{Z}_{p}[X]$.
    '''
    assert has_int_coeffs(f)
    assert has_int_coeffs(g)

    gcd = poly_mod_gcd(f, g, modulus)
    if not is_const_poly(gcd):
        return False
    const_coeff = gcd[0]
    return integer_ring.is_unit(const_coeff, modulus)


def normalize_to_monic(p: np.ndarray[int], modulus: int) -> np.ndarray[int]:
    r'''
    Normalize a polynomial to monic polynomial by multiplying it with a modular inverse of its leading coefficient.
    '''
    assert has_int_coeffs(p)

    leading_term = p[deg(p, modulus)]
    assert integer_ring.is_unit(leading_term, modulus)

    return poly_mod_scalar_mul(p, integer_ring.modular_inverse(leading_term, modulus))
    

def poly_mod_eea(f: np.ndarray[int], g: np.ndarray[int], modulus: int) -> tuple[np.ndarray[int],np.ndarray[int],np.ndarray[int]]:
    r'''
    Extended Euclidean algorithm on polynomial in ring $\mathcal{Z}_{p}[X]$.
    '''
    assert has_int_coeffs(f)
    assert has_int_coeffs(g)

    f0, f1 = f, g
    a0, a1 = np.array([1], dtype=int), np.array([0], dtype=int)
    b0, b1 = np.array([0], dtype=int), np.array([1], dtype=int)
    while not is_zero_poly(f1, modulus):
        q, r = poly_mod_euclidean_division(f0, f1, modulus)
        f0, f1 = f1, r

        a0, a1 = a1, poly_mod_sub(a0, poly_mod_mul(q, a1, modulus), modulus)
        b0, b1 = b1, poly_mod_sub(b0, poly_mod_mul(q, b1, modulus), modulus)
    return f0, a0, b0


def center_lift(coeffs: np.ndarray[int], modulus: int) -> np.ndarray:
    r'''
    Implementation of center lift that is used in NTRU cryptosystem.
    '''
    assert has_int_coeffs(coeffs)
    
    coeffs_set = np.arange(modulus // 2, -modulus // 2, -1)[::-1]
    coeffs_map = { p % modulus : p for p in coeffs_set}
    return np.array([coeffs_map[c] for c in coeffs], dtype=int)