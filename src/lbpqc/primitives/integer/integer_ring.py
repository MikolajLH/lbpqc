import math
import numpy as np
from lbpqc.type_aliases import *
from typing import Tuple

def eea_int(a: int, b: int) -> tuple[int, int, int]:
    r'''
    Function implementing extended Euclidean algorithm on integers.

    It returns the gcd of a and b, as well as the coefficients (x, y) such that:
    $$ ax + by = gcd(a, b) $$
    '''
    if b == 0:
        return a, 1, 0
    
    gcd, x1, y1 = eea_int(b, a % b)
    x = y1
    y = x1 - (a // b) * y1
    return gcd, x, y


def modular_inverse(a: int, modulus: int):
    r'''
    Function finding modular multiplicative inverse $a'$ of $a$
    such that $a' a \equiv 1 \mod(modulus)$

    if the inverse does not exisit - $\gcd(a,modulus) \neq 1 $ raises an exception

    uses extended Euclidean algorithm
    '''
    gcd, a_inv, _ = eea_int(a, modulus)
    if gcd != 1:
        raise ValueError(f"Modular inverse does not exist because gcd({a}, {modulus}) != 1")
    return a_inv % modulus


def gcd_int(a: int, b: int) -> int:
    r'''
    computes gcd of two integers
    '''
    return math.gcd(a,b)


def is_unit(a: int, modulus: int) -> bool:
    r'''
    Checks if element of a Integer Ring is a unit
    '''
    return gcd_int(a, modulus) == 1


def LWR_rounding(a: int, p: int, q: int):
    r'''
    from Zq to Zp
    '''
    return math.floor((p/q) * a) % p
LWR_rounding = np.vectorize(LWR_rounding)


def center_reduce(a: int, modulus: int):
    r'''
    reduce integer a into interval (-modulus/2, modulus/2] ^ Z
    '''
    return (a % modulus) - math.floor((modulus - 1) / 2)
center_reduce = np.vectorize(center_reduce)



def mod_pow(a: int, n: int, modulus: int):
    assert n >= 0
    y, z = 1, a
    while n != 0:
        if n % 2 == 1: y = (y * z) % modulus
        n //= 2
        z = (z * z) % modulus
    return y



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
    return modular_inverse(a, modulus)


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


def center_lift(a: int, modulus: int) -> int:
    return (a % modulus) - math.floor((modulus - 1) / 2)
center_lift = np.vectorize(center_lift)

