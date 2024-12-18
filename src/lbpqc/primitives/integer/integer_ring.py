from typing import Tuple
import math

from lbpqc.type_aliases import *


@np.vectorize
def LWR_rounding(a: int, q: int, p: int) -> ModInt:
    r'''
    from Zq to Zp
    '''
    return math.floor((p/q) * a) % p


@np.vectorize
def mod_reduce(a: int, modulus: int) -> ModInt:
    return a % modulus


@np.vectorize
def center_mod_reduce(a: int, modulus: int, right_inclusive: bool = True) -> CenteredModInt:
    r'''
    reduce integer $a$ into interval (-modulus/2, modulus/2] if right_inclusive is True
    else reduce integer $a$ into interval [-modulus/2, modulus/2) 
    '''
    if right_inclusive:
        return ((a + modulus // 2) % modulus) - modulus // 2
    else:
        return ((a + 1 + modulus // 2) % modulus) - modulus // 2 - 1


def eea(a: int, b: int) -> Tuple[int,int,int]:
    old_s, s = 1, 0
    old_r, r = a, b
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
    
    t = 0 if b == 0 else (old_r - old_s * a) // b
    s = old_s
    gcd = old_r
    return gcd, s, t


def modinv(a: int, modulus: int) -> ModInt:
    gcd, a_inv, _ = eea(a, modulus)
    if gcd != 1:
        raise ValueError(f"Modular inverse of {a} mod {modulus} does not exist gcd is equal to {gcd}")
    
    return a_inv % modulus


def modpow(a: int, r: int, modulus: int) -> ModInt:
    if r < 0:
        return modpow(modinv(a, modulus), -r, modulus)
    
    y, z = 1, a
    while r != 0:
        if r % 2 == 1:
            y = (y * z) % modulus
        r //= 2
        z = (z * z) % modulus
    return y



# class ModIntRing:
#     def __init__(self, modulus: int) -> None:
#         pass

#     def reduce(self, a: int) -> ModInt:
#         pass

#     def center_reduce(self, a: int) -> CenteredModInt:
#         pass

#     def add(self, a, b) -> ModInt:
#         pass

#     def sub(self, a, b) -> ModInt:
#         pass

#     def inv(self, a) -> ModInt:
#         pass