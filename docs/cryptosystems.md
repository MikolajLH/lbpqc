# Cryptosystems submodule

The `cryptosystems` submodule is a collection of **Public Key Cryptosystems** implemented using the `lbpqc.primitives` submodule.

---
Implemented cryptosystems:

* [Congruential pkc](cryptosystems/congruential_pkc.md)
* [Discret knapsack pkc](cryptosystems/discret_knapsack_pkc.md)
* [GGH pkc](cryptosystems/GGH_pkc.md)
* [NTRUEncrypt pkc](cryptosystems/NTRUEncrypt_pkc.md)
* [LWE pkc](cryptosystems/LWE_pkc.md)
* [LWR pkc](cryptosystems/LWR_pkc.md)
* [ringLWE pkc](cryptosystems/ringLWE_pkc.md)
* [ringLWR pkc](cryptosystems/ringLWR_pkc.md)
* [RSA pkc](cryptosystems/RSA_pkc.md)

---

Each cryptosystem is contained in a class with a following API


```python
class Cryptosystem:
    def __init__(self, params):
        '''
        sets up the cryptosytem with params
        '''

    def create_key(self, *params) -> (private_key, public_key):
        '''
        given explicite key parameters creates pair (private_key, public_key)
        '''
        pass

    #TODO
    def generate_key(self) -> (private_key, public_key):
        '''
        generates pair (private_key, public_key) using RNG instance
        '''
        pass

    def encrypt(public_key, plaintext, params) -> ciphertext:
        '''
        Transforms element of Plaintext Space into element of Ciphertext Space, where 
        Plaintext Space and Ciphertext Space are defined for specific cryptosystem
        '''
        pass
    
    def decrypt(private_key, ciphertext) -> plaintext:
        '''
        Transforms element of Ciphertext Space into element of Plaintext Space, where
        Plaintext Space and Ciphertext Space are defined for specific cryptosystem
        '''
        pass

    # TODO
    def plaintext_from_bytes():
        assert(False)

    # TODO
    def bytes_to_plaintext():
        assert(False)

    #TODO
    def encrypt_bytes(public_key, plaintext: bytes, params) -> ciphertext:
        assert(False)

    #TODO
    def decrypt_bytes(private_key, ciphertext) -> plaintext : bytes:
        assert(False)
```