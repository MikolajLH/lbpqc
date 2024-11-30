from typing import Tuple

from lbpqc.type_aliases import *
import lbpqc.primitives.polynomial.poly as poly
from lbpqc.primitives.integer.integer_ring import modinv, center_lift


class ModIntPolyRing:
    def __init__(self, modulus) -> None:
        if modulus <= 1: raise ValueError("Modulus has to be greater than 1")
        self.modulus = modulus

    
    def reduce(self, polynomial: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial): raise TypeError()

        return poly.trim(polynomial % self.modulus)
    

    def is_zero(self, polynomial: VectorInt) -> bool:
        if not is_VectorInt(polynomial): raise TypeError()

        return poly.is_zero_poly(self.reduce(polynomial))
    

    def deg(self, polynomial: VectorInt) -> int:
        if not is_VectorInt(polynomial): raise TypeError()

        return poly.deg(self.reduce(polynomial))


    def add(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

        return self.reduce(poly.add(polynomial_a, polynomial_b))


    def sub(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

        return self.reduce(poly.sub(polynomial_a, polynomial_b))
        

    def mul(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

        return self.reduce(poly.mul(polynomial_a, polynomial_b))
        

    def euclidean_div(self,  polynomial_a: VectorInt, polynomial_b: VectorInt) -> Tuple[VectorMod, VectorMod]:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

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

    def rem(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

        if self.is_zero(polynomial_b): raise ZeroDivisionError("Can't divide by zero polynomial")

        _, r = self.euclidean_div(polynomial_a, polynomial_b)
        return r


    def to_monic(self, polynomial: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial): raise TypeError()

        leading_coeff = polynomial[self.deg(polynomial)]

        return self.reduce(modinv(leading_coeff, self.modulus) * polynomial)


    def gcd(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()


        r0 = self.reduce(polynomial_a)
        r1 = self.reduce(polynomial_b)
        if poly.deg(r1) > poly.deg(r0):
            r0, r1 = r1, r0
        
        while not self.is_zero(r1):
            r0, r1 = r1, self.rem(r0, r1)
        
        return r0

    def coprime(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> bool:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

        return np.all(self.to_monic(self.gcd(polynomial_a, polynomial_b)) == poly.monomial(1, 0))
    

    def eea(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> Tuple[VectorMod, VectorMod, VectorMod]:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

        f0, f1 = self.reduce(polynomial_a), self.reduce(polynomial_b)
        a0, a1 = poly.monomial(1, 0), poly.zero_poly()
        b0, b1 = poly.zero_poly(), poly.monomial(1, 0)

        while not self.is_zero(f1):
            q, r = self.euclidean_div(f0, f1)

            f0, f1 = f1, r

            a0, a1 = a1, self.sub(a0, self.mul(q, a1))
            b0, b1 = b1, self.sub(b0, self.mul(q, b1))
    
        return f0, a0, b0