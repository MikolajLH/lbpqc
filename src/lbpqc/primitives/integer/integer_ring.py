import math


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