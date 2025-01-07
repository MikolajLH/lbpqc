from lbpqc.type_aliases import *
import lbpqc.primitives.polynomial.poly as poly
from lbpqc.primitives.integer.integer_ring import modinv


class ModIntPolyRing:
    r'''This class implements operations over $\mathbb{Z}_p[X]$ polynomial ring, i.e. polynomials which coefficients are reduced modulo $p$.

    It's important to note that polynomials are represented as numpy's arrays with entries corresponding to coefficients in order of increasing powers.  
    Instance of `ModIntPolyRing` class represents the polynomials ring.  
    It's methods implements operations in this rings, but polynomials in their inputs and outputs are still numpy's arrays.

    Attributes:
        modulus (int): modulus $p$ of the $\mathbb{Z}_p[X]$ polynomial ring.
    '''
    @enforce_type_check
    def __init__(self, modulus: int) -> None:
        r'''
        Constructs the ring object for a given **modulus**.

        Args:
            modulus: Ring modulus.
        '''
        if modulus <= 1: raise ValueError("Modulus has to be greater than 1")
        self.modulus = modulus

    
    @enforce_type_check
    def reduce(self, polynomial: VectorInt) -> VectorModInt:
        r'''
        Reduces the given polynomial to it's cannonical equivalence class in the ring, i.e. reduces polynomials' coefficients modulo **modulus**.

        Args:
            polynomial: Polynomial's coefficients array.
        
        Returns:
            Array of coefficients reduced modulo **modulus**.
        '''
        return poly.trim(polynomial % self.modulus)
    

    @enforce_type_check
    def is_zero(self, polynomial: VectorInt) -> bool:
        r'''
        Checks whether given polynomial is zero polynomial in the ring,
        so whether all it's coefficients are divisible by **modulus**.

        Args:
            polynomial: Coefficients array of polynomial to be checked.
        
        Returns:
            `True` if polynomial is zero polynomial in the ring else `False`
        '''
        return poly.is_zero_poly(self.reduce(polynomial))
    

    @enforce_type_check
    def deg(self, polynomial: VectorInt) -> int:
        r'''Calculates degree of the given polynomial.

        Args:
            polynomial: polynomial's coefficients
        
        Returns:
            degree of the polynomial in the ring.
        '''

        return poly.deg(self.reduce(polynomial))


    @enforce_type_check
    def add(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:
        r'''Adds polynomial $a$ to polynomial $b$.

        Args:
            polynomial_a: polynomial's $a$ coefficients.
            polynomial_b: polynomial's $b$ coefficients.
    
        Returns:
            Coefficients array of polynomial $a + b$.
        '''

        return self.reduce(poly.add(polynomial_a, polynomial_b))

    @enforce_type_check
    def sub(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:
        r'''Subtract polynomial $b$ from polynomial $a$.

        Args:
            polynomial_a: polynomial's $a$ coefficients.
            polynomial_b: polynomial's $b$ coefficients.
    
        Returns:
            Coefficients array of polynomial $a - b$.
        '''
        

        return self.reduce(poly.sub(polynomial_a, polynomial_b))
        

    @enforce_type_check
    def mul(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:
        r'''Multiplies polynomial $a$ and polynomial $b$.

        Args:
            polynomial_a: polynomial's $a$ coefficients.
            polynomial_b: polynomial's $b$ coefficients.
    
        Returns:
            Coefficients array of polynomial $a \cdot b$.
        '''

        return self.reduce(poly.mul(polynomial_a, polynomial_b))


    @enforce_type_check
    def euclidean_div(self,  polynomial_a: VectorInt, polynomial_b: VectorInt) -> Tuple[VectorModInt, VectorModInt]:
        r'''Euclidean division (long division) for polynomials in the ring.

        Args:
            polynomial_a: .
            polynomial_b: .

        Returns:
            Quotient and Remainder polynomials.
        '''

        if self.is_zero(polynomial_b): raise ZeroDivisionError("Can't divide by zero polynomial")

        q = poly.zero_poly()
        r = self.reduce(polynomial_a)

        d = self.deg(polynomial_b)
        c = polynomial_b[d]
        while (dr := self.deg(r)) >= d:
            s = poly.monomial(r[dr] * modinv(c, self.modulus), dr - d)
            q = self.add(q, s)
            r = self.sub(r, self.mul(s, polynomial_b))
        
        return q, r


    @enforce_type_check
    def rem(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:
        r'''Computes only remainder of the euclidean divions of two ring's polynomials.

        Args:
            polynomial_a: .
            polynomial_b: .

        Returns:
            Remainder of the euclidean division.
        '''

        if self.is_zero(polynomial_b): raise ZeroDivisionError("Can't divide by zero polynomial")

        _, r = self.euclidean_div(polynomial_a, polynomial_b)
        return r

    @enforce_type_check
    def to_monic(self, polynomial: VectorInt) -> VectorModInt:
        r'''Reduces given polynomial to monic polynomial by multiplying it by scalar equal to modular multiplicative inverse of the highest power coefficient.

        Args:
            polynomial: .

        Returns:
            Polynomial with leading coefficient equal to 1.
        '''
    
        leading_coeff = polynomial[self.deg(polynomial)]

        return self.reduce(modinv(leading_coeff, self.modulus) * polynomial)


    @enforce_type_check
    def gcd(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:
        r'''calculates $\gcd$ of polynomials $a$ and $b$ using euclidean algorithm.  
        It's worth noting that int the polynomial ring, $\gcd$ is not unique.

        Args:
            polynomial_a: Coefficients array of polynomial $a$.
            polynomial_b: Coefficients array of polynomial $b$.

        Returns:
            Polynomial with maximal degree that divides both $a$ and $b$.

        '''
    
        r0 = self.reduce(polynomial_a)
        r1 = self.reduce(polynomial_b)
        if poly.deg(r1) > poly.deg(r0):
            r0, r1 = r1, r0
        
        while not self.is_zero(r1):
            r0, r1 = r1, self.rem(r0, r1)
        
        return r0


    @enforce_type_check
    def coprime(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> bool:
        r'''checks whether two polynomials $a$ and $b$ are coprime in a ring.  
        Checks whether $\gcd(a,b)$ reduced to monic polynomial is equal to constant polynomial $p \equiv 1$.

        Args:
            polynomial_a: Coefficients array of polynomial $a$.
            polynomial_b: Coefficients array of polynomial $b$.

        Returns:
            `True` if $a$ and $b$ are coprime in the ring, `False` otherwise.

        '''
        return np.all(self.to_monic(self.gcd(polynomial_a, polynomial_b)) == poly.monomial(1, 0))
    

    @enforce_type_check
    def eea(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> Tuple[VectorModInt, VectorModInt, VectorModInt]:
        r'''Extended Euclidean algorithm for polynomials, e.i algorithm that calculates coefficients for **Bézout's identity**.

        Args:
            polynomial_a: Coefficients array of polynomial $a$.
            polynomial_b: Coefficients array of polynomial $b$.

        Returns:
            Tuple of polynomials $(d, s, t)$ that satisfies $\gcd(a,b) = d$ and $s \cdot a + t \cdot b = d$.
        
        '''
        
        f0, f1 = self.reduce(polynomial_a), self.reduce(polynomial_b)
        a0, a1 = poly.monomial(1, 0), poly.zero_poly()
        b0, b1 = poly.zero_poly(), poly.monomial(1, 0)

        while not self.is_zero(f1):
            q, r = self.euclidean_div(f0, f1)

            f0, f1 = f1, r

            a0, a1 = a1, self.sub(a0, self.mul(q, a1))
            b0, b1 = b1, self.sub(b0, self.mul(q, b1))
    
        return f0, a0, b0