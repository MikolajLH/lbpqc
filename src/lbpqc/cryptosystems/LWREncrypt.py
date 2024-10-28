import numpy as np
from lbpqc.primitives.polynomial.LWR_polynomial_ring import *



class Cryptosystem:

    def check_params(self, N: int, q: int, p: int, k: int, t: int) -> bool:
        '''
        t < p < q
        '''
        return True

    def __init__(self, N: int, q: int, p: int, k: int, t: int) -> None:
        assert self.check_params(N, q, k)
        self._params = (N, q, p, k, t)

    @property
    def params(self):
        return self._params

    def create_key(self, A: np.ndarray[np.ndarray[np.ndarray[int]]], s: np.ndarray[np.ndarray[int]]):
        N, q, p, k, t = self.params

        b = np.zeros((k,N), dtype=int)
        for i in range(k):
            for j in range(k):
                b[i] = lwr_poly_add(N, q, b[i], lwr_poly_mul(N, q, A[i][j], s[j]))
            b[i] = lwr_poly_round(b[i], q, p)

        private_key = s
        public_key = (A, b)
        return private_key, public_key

    def generate_key(self):
        assert False

    def encrypt(self, public_key: (np.ndarray[np.ndarray[np.ndarray[int]]], np.ndarray[np.ndarray[int]]), plaintext: np.ndarray[int], sp: np.ndarray[np.ndarray[int]]):
        N, q, p, k, t = self.params
        A, b = public_key
        m = plaintext

        bp = np.zeros((k,N), dtype=int)
        for i in range(k):
            for j in range(k):
                bp[i] = lwr_poly_add(N, q, bp[i], lwr_poly_mul(N, q, A[j][i], sp[j]))
            bp[i] = lwr_poly_round(bp[i], q, p)

        v = mipr.poly_mod_scalar_mul(m, (p+1)//2, p)
        for i in range(k):
            v = lwr_poly_add(N, p, v, lwr_poly_mul(N, p, b[i], sp[i]))
        v = lwr_poly_round(v, p, 2*t)

        ciphertext = (bp, v)
        return ciphertext

    def decrypt(self, private_key: np.ndarray[np.ndarray[int]], ciphertext: (np.ndarray[np.ndarray[int]], np.ndarray[int])):
        N, q, p, k, t = self.params
        s = private_key
        bp, v = ciphertext

        mn = mipr.poly_mod_scalar_mul(v, -p//(2*t), p)
        for i in range(k):
            mn = lwr_poly_add(N, p, mn, lwr_poly_mul(N, p, s[i], bp[i]))
        mn = lwr_poly_round(mn, p, 2)

        return mn

    def encrypt_bytes(self):
        assert False

    def decrypt_bytes(self):
        assert False



