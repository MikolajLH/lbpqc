import numpy as np
from lbpqc.primitives.polynomial.LWE_polynomial_ring import *


class Cryptosystem:

    def check_params(self, N: int, q: int, k: int) -> bool:
        '''
        q prime
        '''
        return True

    def __init__(self, N: int, q: int, k: int) -> None:
        assert self.check_params(N, q, k)
        self._params = (N, q, k)

    @property
    def params(self):
        return self._params

    def create_key(self, A: np.ndarray[np.ndarray[np.ndarray[int]]], s: np.ndarray[np.ndarray[int]], e: np.ndarray[np.ndarray[int]]):
        N, q, k = self.params

        t = e
        for i in range(k):
            for j in range(k):
                t[i] = lwa_poly_add(N, q, t[i], lwe_poly_mul(N, q, A[i][j], s[j]))

        private_key = s
        public_key = (A, t)
        return private_key, public_key

    def generate_key(self):
        assert False

    def encrypt(self, public_key: (np.ndarray[np.ndarray[np.ndarray[int]]], np.ndarray[np.ndarray[int]]), plaintext: np.ndarray[int], r: np.ndarray[np.ndarray[int]], e1: np.ndarray[np.ndarray[int]], e2: np.ndarray[int]):
        N, q, k = self.params
        A, t = public_key
        m = plaintext

        u = e1
        for i in range(k):
            for j in range(k):
                u[i] = lwa_poly_add(N, q, u[i], lwe_poly_mul(N, q, A[j][i], r[j]))

        v = lwa_poly_add(N, q, e2, mipr.poly_mod_scalar_mul(m, (q+1)//2, q))
        for i in range(k):
            v = lwa_poly_add(N, q, v, lwe_poly_mul(N, q, t[i], r[i]))

        ciphertext = (u, v)
        return ciphertext

    def decrypt(self, private_key: np.ndarray[np.ndarray[int]], ciphertext: (np.ndarray[np.ndarray[int]], np.ndarray[int])):
        N, q, k = self.params
        s = private_key
        u, v = ciphertext

        mn = mipr.poly_mod_scalar_mul(v, -1, q)
        for i in range(k):
            mn = lwa_poly_add(N, q, mn, lwe_poly_mul(N, q, s[i], u[i]))
        mn = mipr.poly_mod_scalar_mul(mn, -1, q)

        for i in range(N):
            mn[i] = ((mn[i]+(q+1)//4)%q)//((q+1)//2)
        return mn

    def encrypt_bytes(self):
        assert False

    def decrypt_bytes(self):
        assert False