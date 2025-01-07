from typing import Tuple
import math

from lbpqc.type_aliases import *


@np.vectorize
def LWR_rounding(a: int, q: int, p: int) -> ModInt:
    r'''
    **LWR** rounding function, that maps elements from $\mathbb{Z}_q$ to elements from $\mathbb{Z}_p$ for some positive integers $p$ and $q$.

    Args:
        a: Integer to be rounded.
        q: Modulus of domain ring $\mathbb{Z}_q$.
        p: Modulus of codomain ring $\mathbb{Z}_p$.
    
    Returns:
        Integer $c \in [0,p)$.
    '''
    return math.floor((p/q) * a) % p


@np.vectorize
def mod_reduce(a: int, m: int) -> ModInt:
    r'''
    Reduces integer $a$ to it's equivalence class modulo $m$ represented as integer in the interval $[0,m)$.
    
    Args:
        a: integer or numpy array to be reduced.
        m: positive modulus for the congruence relation.
    
    Returns:
        Integer or numpy array with entries from interval $[0, \text{m})$.
    '''
    return a % m


@np.vectorize
def center_mod_reduce(a: int, m: int, right_closed: bool = True) -> CenteredModInt:
    r'''
    Reduces integer $a$ to it's equivalence class modulo $m$ represented as an interval **centered around zero**.
    Depending on the `right_closed` parameter, the interval is either
    $$
    \left(-\frac{m}{2}, \frac{m}{2}\right]
    $$
    or
    $$
    \left[-\frac{m}{2}, \frac{m}{2}\right)
    $$

    Args:
        a: integer or numpy array to be reduced.
        m: positive modulus for the congruence relation.
        right_closed: parameter deciding which side of half-open interval is closed.
    
    Returns:
        Integer or numpy array with entries reduced around zero.
    '''
    if right_closed:
        return ((a + m // 2) % m) - m // 2
    else:
        return ((a + 1 + m // 2) % m) - m // 2 - 1


def eea(a: int, b: int) -> Tuple[int,int,int]:
    r'''
    Implementation of **extended Euclidean algorithm**, i.e. algorithm that for integers $a$ and $b$ computes their $\gcd$ and
    coefficients of **Bézout's identity**, which are integers $s$ and $t$ such that
    $$
    sa + tb = \gcd(a,b)
    $$

    Args:
        a: integer a.
        b: integer b.
    
    Returns:
        Tuple (gcd(a,b), s, t).
    '''
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
    r'''Computes multiplicative modular inverse of integer $a$ for a given modulus.
    If the inverse does not exists, i.e. $\gcd(a, \text{modulus}) \neq 1$ then raises `ValueError`.

    Args:
        a: integer, which inverse we want to calculate.
        modulus: positive modulus.
    
    Returns:
        Integer $r$ from interval $[0, \text{modulus})$ that satisfies $a \cdot r \equiv 1 \mod \text{modulus}$.
    
    Raises:
        ValueError: If $\gcd(a, \text{modulus}) \neq 1$.
    '''
    gcd, a_inv, _ = eea(a, modulus)
    if gcd != 1:
        raise ValueError(f"Modular inverse of {a} mod {modulus} does not exist gcd is equal to {gcd}")
    
    return a_inv % modulus


def modpow(a: int, r: int, modulus: int) -> ModInt:
    r'''Computes
    $$
    a^{r} \mod \text{modulus}
    $$
    using *multiply and halve* powering algorithm for groups.

    Args:
        a: Base.
        r: Exponent.
        modulus: Positive modulus.
    
    Returns:
        Integer $c \in [0, \text{modulus})$ such that $c \equiv a^r \mod \text{modulus}$.

    '''
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