import numpy as np
from lbpqc.primitives.polynomial.NTRU_polynomial_ring import *
from lbpqc.primitives.polynomial import random



class Cryptosystem:

    def check_params(self, N: int, p: int, q: int, d: int) -> bool:
        '''
        N prime
        p prime
        gcd(p,q) = gcd(N, q) = 1
        q > (6d + 1)p
        '''
        return True
        

    def __init__(self, N: int, p: int, q: int, d: int) -> None:
        assert self.check_params(N, p, q, d)
        self._params = (N, p, q, d)

    @property
    def params(self):
        return self._params
    
    def create_key(self, f: np.ndarray[int], g: np.ndarray[int]):
        N, p, q,d = self.params

        Fq = ntru_poly_inverse(N, q, f)
        Fp = ntru_poly_inverse(N, p, f)

        h = ntru_poly_mul(N, q, Fq, g)

        private_key = (f, Fp)
        public_key = h
        return private_key, public_key
    

    def generate_key(self):
        N, p, q, d = self.params

        f = random.sample_ternary_polynomial_coeffs(N, d + 1, d)
        g = random.sample_ternary_polynomial_coeffs(N, d, d)

        return self.create_key(f, g)


    def encrypt(self, public_key: np.ndarray[int], plaintext: np.ndarray[int], noise: np.ndarray[int]):
        h = public_key
        m = plaintext
        r = noise

        N, p, q, d = self.params

        e = mipr.poly_mod_add(ntru_poly_mul(N, q, mipr.poly_mod_scalar_mul(h, p, modulus=q), r), m, q)
        return e

    def decrypt(self, private_key, ciphertext: np.ndarray[int]):
        N, p, q, d = self.params
        f, Fp = private_key
        e = ciphertext

        a = mipr.center_lift(ntru_poly_mul(N, p, Fp, mipr.center_lift(ntru_poly_mul(N, q, f, e), q)), p)
        return mipr.trim_trailing_zeros(a, p)


    def encrypt_bytes(self):
        assert False

    def decrypt_bytes(self):
        assert False