from lbpqc.cryptosystems.NTRUEncrypt_pkc import Cryptosystem
import numpy as np 


N, p, q, d = 7, 3, 41, 2
f = np.array([-1, 0, 1, 1, -1, 0, 1], dtype=int)
g = np.array([0, -1, -1, 0, 1, 0, 1], dtype=int)

Fq = np.array([37, 2, 40, 21, 31, 26, 8], dtype=int)
Fp = np.array([1, 1, 1, 1, 0, 2, 1], dtype=int)

h = np.array([30, 26, 8, 38, 2, 40, 20], dtype=int)

m = np.array([ 1,-1, 1, 1, 0,-1], dtype=int)
r = np.array([-1, 1, 0, 0, 0,-1, 1], dtype=int)

e = np.array([25, 3, 40, 2, 4, 19, 31], dtype=int)

cs = Cryptosystem(N,p,q,d)


def test_create_key():
    priv, publ = cs.create_key(f, g)
    assert (priv[0] == f).all()
    assert (priv[1] == Fp).all()
    assert (publ == h).all()


def test_encrypt():
    cipher = cs.encrypt(h, m, r)
    assert (cipher == e).all()

def test_decrypt():
    decr = cs.decrypt((f, Fp), e)
    assert (decr == m).all()