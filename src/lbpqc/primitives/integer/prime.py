from . import integer_ring
import random


# All primes less than 256
SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251)


def fermat_primality_test(p: int, s: int, int_gen):
    r'''
    Fermat's primality test
    '''
    if p in SMALL_PRIMES:
        return True
    
    for _ in range(s):
        a = int_gen(2, p - 2)
        if integer_ring.modpow(a, p - 1, p) == 1:
            return False
    return True


def miller_rabin_primality_test(p: int, s: int, int_gen):
    '''
    Miller-Rabin primality test
    '''
    if p in SMALL_PRIMES:
        return True
    
    if p == 1 or p % 2 == 0 or p % 3 == 0 or p % 5 == 0:
        return False
    u = 0
    r = p - 1
    while r % 2 == 0:
        u += 1
        r //= 2

    for _ in range(s):
        a = int_gen(2, p - 2)
        z = integer_ring.modpow(a, r, p)
        if z != 1 and z != p - 1:
            for _ in range(u - 1):
                z = (z * z) % p
                if z == 1:
                    return False
            if z != p - 1:
                return False
        return True
    

def is_prime(p: int) -> bool:

    if p in SMALL_PRIMES:
        return True
    
    return miller_rabin_primality_test(p, 20, lambda a, b: random.randint(a, b - 1))