import itertools
from lbpqc.primitives.integer import integer_ring

class Cryptosystem:
    def __init__(self) -> None:
        pass

    def is_superincreasing(self,rs: list[int]):
        return all((j >= 2 * i for i,j in itertools.pairwise(rs)))
    
    def assert_key_parameters(self, rs: list[int], A: int, B: int):
        assert self.is_superincreasing(rs)
        assert B > 2 * rs[-1]
        assert integer_ring.gcd_int(A, B) == 1
    
    def create_key(self, rs: list[int], A: int, B: int):
        self.assert_key_parameters(rs, A, B)

        Ms = [(A * ri) % B for ri in rs]

        private_key = (rs, A, B)
        public_key = Ms

        return private_key, public_key

    def generate_key(self):
        assert False

    def encrypt(self, Ms: list[int], xs: list[bool]) -> int:
        assert len(xs) <= len(Ms)

        S = sum((x * M for x, M in zip(xs, Ms) if x))
        return S

    def decrypt(self, rs: list[int], A: int, B: int, S: int):
        Sprim = (integer_ring.modular_inverse(A, B) * S) % B
        xs = [0 for _ in rs]
        for i,r in enumerate(rs[::-1]):
            if r <= Sprim:
                xs[i] = 1
                Sprim = Sprim - r

        return xs[::-1]

    def encrypt_bytes(self, public_key, plaintext: bytes) -> int:
        assert False

    def decrypt_bytes(self, private_key, ciphertext: int) -> bytes:
        assert False