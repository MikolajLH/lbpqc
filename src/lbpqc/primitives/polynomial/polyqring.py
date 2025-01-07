from lbpqc.type_aliases import *

from lbpqc.primitives.integer import integer_ring
from lbpqc.primitives.polynomial import poly, modpoly


class PolyQuotientRing:
    r'''This class implements operations over polynomial quotient ring
    $$
    \frac{\mathbb{Z}_p[X]}{q(X)}
    $$
    for given positive integer $p$ and polynomial $q$.

    It's important to note that polynomials are represented as numpy's arrays with entries corresponding to coefficients in order of increasing powers.  
    Instance of `PolyQuotientRing` class represents the polynomials ring.  
    It's methods implements operations in this rings, but polynomials in their inputs and outputs are still numpy's arrays.

    Attributes:
        poly_modulus (VectorInt): .
        int_modulus (int): .
        Zm (ModIntPolyRing): Object representing $\mathbb{Z}_p$ ring.
    '''
    @enforce_type_check
    def __init__(self, poly_modulus: VectorInt, int_modulus: int) -> None:
        r'''Constructs the ring object for a given polynomial modulus and integer modulus.

        Args:
            poly_modulus: Coefficients array of polynomial modulus fot the ring. The $q(X)$ in $\frac{\mathbb{Z}_p[X]}{q(X)}$.
            int_modulus: Integer modulus for the ring. The $p$ in $\frac{\mathbb{Z}_p[X]}{q(X)}$.
        '''
        self.poly_modulus = poly_modulus
        self.int_modulus = int_modulus
        self.Zm = modpoly.ModIntPolyRing(int_modulus)

    
    @property
    def quotient(self):
        return self.poly_modulus
    

    @enforce_type_check
    def reduce(self, polynomial: VectorInt) -> VectorModInt:
        r'''Reduces the given polynomial $u$ to it's cannonical equivalence class in the ring,
        i.e. takes, the remainder of division $\frac{u}{q}$, where $q$ is polynomial modulus for the ring.  
        The division is performed in $\mathbb{Z}_p[X]$ ring.

        Args:
            polynomial: Polynomial's coefficients array.
        
        Returns:
            Coefficients array of polynomial that is the remainder of the Euclidean division.
        '''
        return self.Zm.rem(polynomial, self.poly_modulus)
        

    @enforce_type_check
    def add(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:
        r'''Adds polynomial $a$ to polynomial $b$.

        Args:
            polynomial_a: polynomial's $a$ coefficients.
            polynomial_b: polynomial's $b$ coefficients.
    
        Returns:
            Coefficients array of polynomial $a + b$.
        '''

        return self.reduce(self.Zm.add(polynomial_a, polynomial_b))
    

    @enforce_type_check
    def sub(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:
        r'''Subtract polynomial $b$ from polynomial $a$.

        Args:
            polynomial_a: polynomial's $a$ coefficients.
            polynomial_b: polynomial's $b$ coefficients.
    
        Returns:
            Coefficients array of polynomial $a - b$.
        '''

        return self.reduce(self.Zm.sub(polynomial_a, polynomial_b))
    
    @enforce_type_check
    def mul(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:
        r'''Multiplies polynomial $a$ and polynomial $b$.

        Args:
            polynomial_a: polynomial's $a$ coefficients.
            polynomial_b: polynomial's $b$ coefficients.
    
        Returns:
            Coefficients array of polynomial $a \cdot b$.
        '''
        
        return self.reduce(self.Zm.mul(polynomial_a, polynomial_b))
        
    @enforce_type_check
    def inv(self, polynomial: VectorInt) -> VectorModInt:
        r'''For a given polynomial $v$, calculates it's multiplicative inverse in a ring.

        Args:
            polynomial: Polynomial's coefficients array.
        
        Returns:
            Coefficients array of polynomial $u$, such that $u \cdot v \equiv 1$.

        Raises:
            ValueError: When given polynomial is not an unit in the ring (it's not coprime with ring's polynomial modulus).

        '''
        
        if not self.Zm.coprime(polynomial, self.poly_modulus): raise ValueError("Inverse does not exists")

        gcd, u, _ = self.Zm.eea(polynomial, self.poly_modulus)

        c = integer_ring.modinv(gcd, self.int_modulus)

        return self.reduce(u * c)


def construct_ring(p: str, N: int, q: int) -> PolyQuotientRing|None:
    r'''Function for constructing commonly used quotient rings.

     Possible optinos:  
        - "-" || "X^N - 1": represents $X^N - 1$ polynomial family.  
        - "+" || "X^N + 1": represents $X^N + 1$ polynomial family.  
        - "prime" || "X^N - x - 1": represents $X^N - X - 1$ polynomial family.


    Args:
        p: string representation/symbol of ring's polynomial modulus.
        N: Degree of the polynomial modulus of the ring.
        q: Integer modulus of the ring.

    Returns:
        None if the parameters were invalid else ring object.
    '''
    g = poly.zero_poly(N)
    match p:
        case "-" | "X^N - 1":
            g[[0, N]] =-1, 1
            pass
        case "+" | "X^N + 1":
            g[[0, N]] = 1, 1
            pass
        case "prime" | "X^N - x - 1":
            g[[0, 1, N]] = -1, -1, 1
        case _:
            return None
        
    return PolyQuotientRing(g, q)