from lbpqc.cryptosystems.congruential_pkc import Cryptosystem

q = 122430513841
f = 231231
g = 195698
h = 39245579300

m = 123456
r = 101010
e = 18357558717
d = 123456

cs = Cryptosystem(q)

def test_create_key():
    priv, publ = cs.create_key(f, g)
    assert priv == (f,g)
    assert publ == h


def test_encrypt():
    cipher = cs.encrypt(h, m, r)
    assert cipher == e


def test_decrypt():
    decr = cs.decrypt(f, g, e)
    assert decr == d