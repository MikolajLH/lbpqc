from lbpqc.cryptosystems.LWEEncrypt_pkc import Cryptosystem
import numpy as np


N, q, k = 4, 17, 2

A = np.array([np.array([np.array([11,16,16,6]),np.array([3,6,4,9])]),np.array([np.array([1,10,3,5]),np.array([15,9,1,6])])])
s = np.array([np.array([0,1,-1,-1]),np.array([0,-1,0,-1])])
e = np.array([np.array([0,0,1,0]),np.array([0,-1,1,0])])

t = np.array([np.array([7,0,15,16]),np.array([6,11,12,10])])

m = np.array([1,1,0,1])
r = np.array([np.array([0,0,1,-1]),np.array([-1,0,1,1])])
e1 = np.array([np.array([0,1,1,0]),np.array([0,0,1,0])])
e2 = np.array([0,0,-1,-1])

u = np.array([np.array([3,10,11,11]),np.array([11,13,4,4])])
v = np.array([16,9,6,8])


cs = Cryptosystem(N,q,k)



def test_create_key():
    priv, publ = cs.create_key(A,s,e)
    assert (priv == s).all()
    assert (publ[0] == A).all()
    assert (publ[1] == t).all()


def test_encrypt():
    cipher = cs.encrypt((A,t), m, r, e1, e2)
    assert (cipher[0] == u).all()
    assert (cipher[1] == v).all()

def test_decrypt():
    decr = cs.decrypt(s, (u,v))
    assert (decr == m).all()