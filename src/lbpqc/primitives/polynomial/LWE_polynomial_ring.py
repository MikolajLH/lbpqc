import numpy as np
import lbpqc.primitives.polynomial.modulo_integer_polynomial_ring as mipr



def reduce_coeffs_by_lwe_modulus_polynomial(N: int, coeffs: np.ndarray[int]):
    r'''
    Reducing powers of polynomial modulo N.
    '''
    assert coeffs.dtype == int
    v = np.pad(coeffs, (0, N - (len(coeffs) % N)))
    # Polynomial is already reduced by (X^N + 1)
    if len(v) == N:
        return v
    m = len(v) // N
    z=[1,-1] * (m // 2)
    if m % 2 == 1:
        z.append(1)
    return ((v.reshape((m, N)))*np.array(z).reshape((m,1))).sum(axis=0)


def lwe_poly_mul(N: int, modulus: int, a_coeffs: np.ndarray[int], b_coeffs: np.ndarray[int]) -> np.ndarray[int]:
    r'''
    Multiplication of two polynomials in in a ring $\frac{\mathcal{Z}_{p}}{{X^N + 1}}\left[X\right]$.
    '''
    return mipr.reduce_coeffs_modulo(
        reduce_coeffs_by_lwe_modulus_polynomial(N, mipr.poly_mod_mul(a_coeffs, b_coeffs, modulus)), modulus)


def lwa_poly_add(N: int, modulus: int, a_coeffs: np.ndarray[int], b_coeffs: np.ndarray[int]) -> np.ndarray[int]:
    r'''
    Addition of two polynomials in a ring $\frac{\mathcal{Z}_{p}}{{X^N + 1}}\left[X\right]$.
    '''
    return mipr.reduce_coeffs_modulo(
        reduce_coeffs_by_lwe_modulus_polynomial(N, mipr.poly_mod_add(a_coeffs, b_coeffs, modulus)), modulus)
