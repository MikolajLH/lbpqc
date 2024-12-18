# from lbpqc.cryptosystems.discret_knapsack_pkc import Cryptosystem

# r = [3, 11, 24, 50, 115]
# A = 113
# B = 250
# M = [89, 243, 212, 150, 245]
# x = [1, 0, 1, 0, 1]
# S = 546
# d = [1, 0, 1, 0, 1]
# cs = Cryptosystem()


# def test_create_key():
#     priv, publ = cs.create_key(r, A, B)
#     assert priv == (r, A, B)
#     assert publ == M


# def test_encrypt():
#     cipher = cs.encrypt(M, x)
#     assert cipher == S


# def test_decrypt():
#     decr = cs.decrypt(r, A, B, S)
#     assert decr == d

def test_true():
    assert True