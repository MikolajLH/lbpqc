from lbpqc.type_aliases import *

import lbpqc.primitives.polynomial.modpoly as modpoly
import lbpqc.primitives.polynomial.poly as poly



class PolyQuotientRing:
    def __init__(self, poly_modulus: VectorInt, int_modulus: int) -> None:
        self.poly_modulus = poly_modulus
        self.int_modulus = int_modulus
        self.Zm = modpoly.ModIntPolyRing(int_modulus)
    

    def reduce(self, polynomial: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial): raise TypeError()
        return self.Zm.rem(polynomial, self.poly_modulus)
        

    def add(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

        return self.reduce(self.Zm.add(polynomial_a, polynomial_b))
    

    def sub(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

        return self.reduce(self.Zm.sub(polynomial_a, polynomial_b))
    

    def mul(self, polynomial_a: VectorInt, polynomial_b: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial_a): raise TypeError()
        if not is_VectorInt(polynomial_b): raise TypeError()

        return self.reduce(self.Zm.mul(polynomial_a, polynomial_b))
        

    def inv(self, polynomial: VectorInt) -> VectorMod:
        if not is_VectorInt(polynomial): raise TypeError()

        if not self.Zm.coprime(polynomial, self.poly_modulus): raise ValueError("Inverse does not exists")

        _, u, _ = self.Zm.eea(polynomial, self.poly_modulus)
        return self.reduce(self.Zm.to_monic(u))
    


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