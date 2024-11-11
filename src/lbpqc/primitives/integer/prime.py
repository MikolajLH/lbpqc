from lbpqc.primitives.integer.integer_ring import mod_pow
import numpy as np



def fermat_primality_test(p: int, s: int):
    if p in [2, 3, 5, 7]:
        return True
    for _ in range(s):
        a = np.random.randint(2, p - 2)
        if mod_pow(a, p - 1, p) == 1:
            return False
    return True


def miller_rabin_primality_test(p: int, s: int):
    if p in [2,3,5,7,11,13]:
        return True
    if p == 1 or p % 2 == 0 or p % 3 == 0 or p % 5 == 0:
        return False
    u = 0
    r = p - 1
    while r % 2 == 0:
        u += 1
        r //= 2

    for _ in range(s):
        a = np.random.randint(2, p - 2)
        z = mod_pow(a, r, p)
        if z != 1 and z != p - 1:
            for _ in range(u - 1):
                z = (z * z) % p
                if z == 1:
                    return False
            if z != p - 1:
                return False
        return True