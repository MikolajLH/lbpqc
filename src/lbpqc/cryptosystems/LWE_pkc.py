from lbpqc.type_aliases import *
from lbpqc.primitives.rng import RNG

class Cryptosystem:
    r'''
    The LWE assumption, with parameters $n, m, q \in \mathbb{N}$ and a "small" error distribution $\chi$ over $\mathbb{Z}$,
    states that for uniformly random
    $$
    \begin{aligned}
    A &~ \mathbb{Z}_{q}^{m \times n}  
    s &~ \mathbb{Z}_{q}^{n}  
    u &~ \mathbb{Z}_{q}^{m}
    \end{aligned}
    $$
    and an error vector $ e ~ \chi^{m} $
    $$
    (A, A \cdot s + e) \text{is computationally indistinguishable from} (A, u)
    $$

    * Learning with rounding revisted; page 1; Introductions; *
    
    Let $q \geq 2$ be an integral modulus.  
    Let $s \in \mathbb{Z}_{q}^{n}$ be an n-dimensional vector.
    Let $\chi$ be a probability distribution over $\mathbb{Z}_{q}$.

    Define $A_{S,\chi}$ as the probability distribution on $\mathbb{Z}_{q}^{n} \times \mathbb{Z}_{q}$,  
    where the samples of this distribution are obtained by the following procedure:
    $$
    \begin{aligned}
    \text{1.} & \text{Take} a \in \mathbb{Z}_{q}^{n} \text{uniformly at random}. \newline
    \text{2.} & \text{Take} e \in \mathbb{Z}_{q} \text{according to the distribution} \chi. \newline
    \text{3.} & \text{Return the tuple} (a, \langle a \; , \; s \rangle + e) \mod q.
    \end{aligned}

    The LWE problem: 
    Given arbitrary number of samples from $A_{s,\chi}$ find $s$.

    * Lattice-based cryptography; page 63; *
    $$

    '''
    def __init__(self, n: int, m: int, q: int, alpha: float, seed: int) -> None:
        self._params = (n, m, q, alpha)
        self._rng = RNG(seed)

    @property
    def params(self):
        return self._params
    
    @property
    def rng(self):
        return self._rng

    
    def create_key(self):
        n, m, q, alpha = self.params
        
        s = self.rng.sample_uniform_Zq(q, n)
        e = self.rng.sample_naive_discrete_gaussian(q, alpha, m)
        A = self.rng.sample_uniform_Zq(q, (n, m))
        b = s @ A + e

        private_key = s
        public_key = (A, b)
        
        return private_key, public_key
    
    
    def generate_key(self):
        return self.create_key()
    

    def encrypt(self, public_key, plaintext: 0|1):
        n, m, q, alpha = self.params
        A, b = public_key
        bit = plaintext

        S = self.rng.random_Zq_subset(q)
        ciphertext = (A.T[S].sum(axis=0), (q//2) * bit + b[S].sum())
        return ciphertext

        # ciphertext = []
        # for bit in plaintext:
        #     S = self.rng.random_Zq_subset(q)
        #     ciphertext.append((A.T[S].sum(axis=0), (q//2) * bit + b[S].sum()))
        #     # if bit == 0:
        #     #     ciphertext.append((A.T[S].sum(axis=0), b[S].sum()))
        #     # else:
        #     #     ciphertext.append((A.T[S].sum(axis=0), q//2 + b[S].sum()))

        # return ciphertext


    
    def decrypt(self, private_key, ciphertext) -> 0|1:
        n, m, q, alpha = self.params
        s = private_key
        a, b = ciphertext
        x = (b - np.dot(s, a)) % q
        if min(x, abs(x - q)) < abs(x - q/2):
            return 0
        else:
            return 1
        

    def encrypt_bytes(self, public_key, plaintext):
        assert False

    def decrypt_bytes(self, private_key, ciphertext):
        assert False