from lbpqc.primitives.lattice import full_rank_lattice as lattice
import numpy as np



class Cryptosystem:
    def __init__(self) -> None:
        self._public_key = None
    
    @property
    def public_key(self):
        return self._public_key

    def create_key(self, good_basis: np.ndarray, U: np.ndarray):
        bad_basis = U @ good_basis

        private_key = good_basis
        self._public_key = bad_basis

        return private_key, self.public_key
    
    def generate_key(self):
        assert False

    def encrypt(self, bad_basis: np.ndarray, plaintext: np.ndarray[int], perturbation_vector: np.ndarray) -> np.ndarray:
        ciphertext = (plaintext @ bad_basis) + perturbation_vector
        return ciphertext

    def decrypt(self, good_basis: np.ndarray, bad_basis: np.ndarray, ciphertext: np.ndarray) -> np.ndarray[int]:
        '''
        z = v @ G
        z = x @ B
        v @ G = x @ B
        v @ (G @ B) = x
        '''
        lattice_vector_in_good_basis = lattice.babai_closest_vector_algorithm(ciphertext, good_basis)
        lattice_vector_in_bad_basis = np.rint(lattice_vector_in_good_basis @ ( good_basis @ np.linalg.inv(bad_basis))).astype(int)
        return lattice_vector_in_bad_basis

    def encrypt_bytes(self, public_key, plaintext: bytes):
        assert False

    def decrypt_bytes(self, private_key, ciphertext) -> bytes:
        assert False


