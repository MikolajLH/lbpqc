from lbpqc.type_aliases import *
from lbpqc.primitives.random import RNG
from lbpqc.primitives.integer.integer_ring import LWR_rounding, center_reduce

class Cryptosystem():
    r'''
    LWR-Based Fully Homomorphic Encryption, Revisited; page 5; 3. building_blocks
    '''
    def __init__(self, n: int, m: int, p: int, q: int, seed: int) -> None:
        self._rng = RNG(seed)
        assert q % p == 0
        self._params = (n, m, p, q)

    @property
    def params(self):
        return self._params

    @property
    def rng(self):
        return self._rng
    
    def create_key(self):
        n, m, p, q = self.params
        sprim =  self.rng.sample_uniform_Zq(2, n)
        A = center_reduce(self.rng.sample_uniform_Zq(q, (n, m)), q)
        b = center_reduce(LWR_rounding(sprim @ A, p, q), p)
        A = np.block([[A], [b]])
        s = np.block([-p/q * sprim, 1])

        public_key = A
        private_key = s

        return private_key, public_key
    
    def encrypt(self, public_key, plaintext: 0|1):
        n, m, p, q = self.params
        A = public_key
        u = plaintext
        r = self.rng.sample_uniform_Zq(2, m)
        c = r @ A
        c[:-1] = center_reduce(c[:-1], q)
        c[-1] = center_reduce(c[-1], p)

        z = np.block([np.zeros(n, dtype=int), u * round(p/2)])
        c = c + z

        c[:-1] = center_reduce(c[:-1], q)
        c[-1] = center_reduce(c[-1], p)

        return c
    
    def decrypt(self, private_key, ciphertext):
        n, m, p, q = self.params
        s = private_key
        c = ciphertext
        u = round((1 / round(p/2)) * np.dot(s,c)) % 2
        
        plaintext = u
        return plaintext
