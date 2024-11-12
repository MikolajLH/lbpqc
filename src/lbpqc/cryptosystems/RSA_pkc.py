from lbpqc.primitives.integer.integer_ring import mod_pow, gcd_int, modular_inverse


class Cryptosystem:
    def __init__(self) -> None:
        pass

    def create_key(self, p: int, q: int, e: int):
        # assert p is prime
        # assert q is prime
        phi = (p - 1) * (q - 1)
        assert 1 <= e <= phi - 1
        assert gcd_int(e, phi) == 1

        n = p * q
        d = modular_inverse(e, phi)

        pk = (n, e)
        sk = d

        return sk, pk
    
    def encrypt(self, public_key, plaintext: int):
        n, e = public_key
        x = plaintext
        y = mod_pow(x, e, n)

        ciphertext = y
        return ciphertext
        

    def decrypt(self, private_key, public_key, ciphertext):
        d = private_key
        n,e = public_key
        y = ciphertext
        x = mod_pow(y, d, n)
        decrypted = x
        return decrypted