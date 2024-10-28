import numpy as np
from lbpqc.primitives.lattice.reduction import *
from lbpqc.primitives.lattice.full_rank_lattice import *
from lbpqc.type_aliases import *



class RNG:
    def __init__(self, seed) -> None:
        r'''
        Initialize numpy rng with seed value.  
        Use `secrets.randbits(128)` for more cryptographicly secure rng.
        '''
        self._rng = np.random.default_rng(seed)
    
    @property
    def rng(self):
        r'''
        Access the underlying `numpy` rng object.
        '''
        return self._rng
    
    
    def sample_naive_discrete_gaussian(self, q: int, alpha: float, size = None) -> int | VectorInt | MatrixInt:
        r'''
        Sample from Discrete Gaussian Distribution defined as
        normal distribution with **mean** zero and **standard deviation** $\alpha q$ that is **rounded to the nearest integer**.  
        $\alpha > 0$ is generally taken such that $\alpha^{-1}$ is a polynomial in $n$, that is $\alpha \approx \frac{1}{n^{c}}$
        for some constant $c$.

        *Lattice-based cryptography; page 64; Definition 3.5 (Learning With Errors (LWE));*
        '''
        return np.rint(self.rng.normal(0, (q * alpha) / (2 * np.pi), size)).astype(int)
    
    
    def sample_uniform_Zq(self, q: int, size : None | int | tuple[int,int] = None) -> ModInt | VectorMod | MatrixMod:
        r'''
        Sample uniformly from $\mathbb{Z}_{q}$ ring.  
        If size is None, returns single element.  
        If size is an int, returns vector (1 dim np.ndarray) with given size.  
        If size is a tuple, returns matrix (2 dim np.ndarray) with given shape.
        '''
        return self.rng.integers(0, q, size)
    
    def random_Zq_subset(self, q: int):
        subset_size = self.rng.integers(0, q)
        return self.rng.choice(q, subset_size, replace=False)