from lbpqc.type_aliases import *

from lbpqc.primitives.integer import integer_ring
from lbpqc.primitives.polynomial import poly, modpoly


class PolyQuotientRing:
    @enforce_type_check
    def __init__(self, poly_modulus: VectorInt, int_modulus: int) -> None:
        self.poly_modulus = poly_modulus
        self.int_modulus = int_modulus
        self.Zm = modpoly.ModIntPolyRing(int_modulus)

    
    @property
    def quotient(self):
        return self.poly_modulus
    

    @enforce_type_check
    def reduce(self, polynomial: VectorInt) -> VectorModInt:
        return self.Zm.rem(polynomial, self.poly_modulus)
        

    @enforce_type_check
    def add(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:

        return self.reduce(self.Zm.add(polynomial_a, polynomial_b))
    

    @enforce_type_check
    def sub(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:

        return self.reduce(self.Zm.sub(polynomial_a, polynomial_b))
    
    @enforce_type_check
    def mul(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorModInt:
        
        return self.reduce(self.Zm.mul(polynomial_a, polynomial_b))
        
    @enforce_type_check
    def inv(self, polynomial: VectorInt) -> VectorModInt:
        
        if not self.Zm.coprime(polynomial, self.poly_modulus): raise ValueError("Inverse does not exists")

        gcd, u, _ = self.Zm.eea(polynomial, self.poly_modulus)

        c = integer_ring.modinv(gcd, self.int_modulus)

        return self.reduce(u * c)


def from_ideal(p: str, N: int, q: int) -> PolyQuotientRing|None:
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