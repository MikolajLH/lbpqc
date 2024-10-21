import numpy as np
import lbpqc.primitives.polynomial.modulo_integer_polynomial_ring as mipr


def reduce_coeffs_by_ntru_modulus_polynomial(N: int, coeffs: np.ndarray[int]):
    r'''
    Reducing powers of polynomial modulo N.
    '''
    assert coeffs.dtype == int
    v = np.pad(coeffs, (0, N - (len(coeffs) % N)))
    # Polynomial is already reduced by (X^N - 1)
    if len(v) == N:
        return v
    m = len(v) // N
    return v.reshape((m, N)).sum(axis=0)


def ntru_poly_mul(N: int, modulus: int, a_coeffs: np.ndarray[int], b_coeffs: np.ndarray[int]) -> np.ndarray[int]:
    r'''
    Multiplication of two polynomials in in a ring $\frac{\mathcal{Z}_{p}}{{X^N - 1}}\left[X\right]$.
    '''
    return mipr.reduce_coeffs_modulo(reduce_coeffs_by_ntru_modulus_polynomial(N, mipr.poly_mod_mul(a_coeffs, b_coeffs, modulus)), modulus)


def ntru_poly_add(N: int, modulus: int, a_coeffs: np.ndarray[int], b_coeffs: np.ndarray[int]) -> np.ndarray[int]:
    r'''
    Addition of two polynomials in a ring $\frac{\mathcal{Z}_{p}}{{X^N - 1}}\left[X\right]$.
    '''
    return mipr.reduce_coeffs_modulo(reduce_coeffs_by_ntru_modulus_polynomial(N, mipr.poly_mod_add(a_coeffs, b_coeffs, modulus)), modulus)
    

def ntru_poly_inverse(N: int, modulus: int, coeffs: np.ndarray[int]) -> np.ndarray[int]:
    r'''
    Calculate multiplicative inverse of a polynomial in a ring $\frac{\mathcal{Z}_{p}}{{X^N - 1}}\left[X\right]$.
    '''
    mpc = mipr.poly_NTRU_modulus(N)
    assert mipr.poly_mod_coprime(coeffs, mpc, modulus)

    gcd, u, _ = mipr.poly_mod_eea(coeffs, mpc, modulus)
    gcd_inv = mipr.integer_ring.modular_inverse(gcd, modulus)

    return mipr.reduce_coeffs_modulo(reduce_coeffs_by_ntru_modulus_polynomial(N,mipr.poly_mod_scalar_mul(u, gcd_inv, modulus)), modulus)