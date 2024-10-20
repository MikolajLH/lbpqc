from lbpqc.primitives.integer import integer_ring

class Cryptosystem:
    def __init__(self, modulus: int):
        self._modulus = modulus
        self.random_element = 101010

    @property
    def modulus(self):
        return self._modulus
    
    def assert_key_parameters(self, f: int, g: int) -> bool:
        q = self.modulus
        lb = (q/4) ** 0.5
        ub = (q/2) ** 0.5
        assert 0 < f < ub
        assert lb < g < ub
        assert integer_ring.gcd_int(f, q) == 1
        assert integer_ring.gcd_int(f, g) == 1
    
    def check_plaintext(self, m: int) -> bool:
        q = self.modulus
        return 0 < m < (q/4) ** 0.5
    
    def check_random_element(self, r: int) -> bool:
        q = self.modulus
        return 0 < r < (q / 2) ** (0.5)
    
    def create_key(self, secret_f: int, secret_g: int):
        self.assert_key_parameters(secret_f, secret_g)

        h = (integer_ring.modular_inverse(secret_f, self.modulus) * secret_g) % self.modulus
        public_key = h
        private_key = (secret_f, secret_g)

        return private_key, public_key
        
    def generate_key(self):
        assert False

    def generate_random_element(self) -> int:
        assert False

    def encrypt(self, h, m, r):
        assert self.check_plaintext(m)
        assert self.check_random_element(r)
        
        e = (r * h + m) % self.modulus
        return e
    
    def decrypt(self, f, g, e):
        a = (f * e) % self.modulus
        b = (integer_ring.modular_inverse(f, g) * a) % g
        return b

    def encrypt_bytes(self, public_key, plaintext: bytes) -> int:
        assert False


    def decrypt_bytes(self, private_key, ciphertext: int) -> bytes:
        assert False