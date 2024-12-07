import numpy as np
from lbpqc.primitives.integer import integer_ring
from lbpqc.type_aliases import *


def elementary_row_swap(n: int, i: int, j: int):
    E = np.identity(n, int)
    E[[i,j]] = E[[j, i]]
    return E

def elementary_row_mul(n: int, i: int, s: int):
    E = np.identity(n, int)
    E[i, i] = s
    return E


def elementary_row_add(n: int, i: int, j: int, s: int):
    E = np.identity(n, int)
    E[i, j] = s
    return E


def HNF(A: MatrixInt) -> Tuple[MatrixInt, SquareMatrixInt, int]:
    r'''
    Computes row-style Hermite Normal Form of a integer matrix A.
    '''
    H = A.copy()
    m, n = H.shape
    p = min(m,n)
    k, j = 0, 0

    U = np.identity(m, dtype=int)
    detU = 1


    while j != p:
        # Choose pivot
        col = H[k:, j]
        non_zero = col[col != 0]
        if len(non_zero) == 0:
            j += 1
            k += 1
            continue
        min_val = np.min(np.abs(non_zero))
        i0 = np.where(np.abs(col) == min_val)[0][0] + k
        if i0 > k:
            H[[k, i0]] = H[[i0, k]]
            detU *= -1
            U = elementary_row_swap(m, k, i0) @ U

        if H[k,j] < 0:
            H[k] = -H[k]
            detU *= -1
            U = elementary_row_mul(m, k, -1) @ U

        # Reduce Rows
        b = H[k,j]
        for i in range(k+1, m):
            q = round(H[i,j] / b)
            H[i] -= q * H[k]
            U = elementary_row_add(m, i, k, -q) @ U

        # Check if column is done
        if np.all(H[k+1:, j] == 0):
            j += 1
            k += 1
            
    # Final reductions
    k = 0
    for j in range(p):
        if H[k,j] < 0:
            H[k] = -H[k]
            U = elementary_row_mul(m, k, -1) @ U
            detU *= -1

        b = H[k,j]
        if b == 0: continue
        for i in range(k):
            q = H[i,j] // b
            H[i] -= q * H[k]
            U = elementary_row_add(m, i, k, -q) @ U

        k += 1
        
    return H, U, detU


def nullity(A: MatrixInt) -> int:
    H, U, _ = HNF(A)
    r = 0
    for row in H[::-1]:
        if np.all(row == 0):
            r += 1
        else:
            break
    return r


def left_kernel(A: MatrixInt) -> MatrixInt|None:
    H, U, _ = HNF(A)
    r = 0
    for row in H[::-1]:
        if np.all(row == 0):
            r += 1
        else:
            break
    if r == 0:
        return None
    return U[-r::]


def det(A: SquareMatrixInt) -> int:
    H, U, detU = HNF(A)
    return np.prod(np.diagonal(H)) * detU


def minor(A: SquareMatrixInt, i: int, j: int) -> int:
    return det(np.delete(np.delete(A, i, axis=0), j, axis=1))


def cofactor(A: SquareMatrixInt, i: int, j: int) -> int:
    return minor(A, i, j) * ((-1) ** (i + 1 + j + 1))


def cofactor_matrix(A: SquareMatrixInt) -> SquareMatrixInt:
    n = A.shape[0]
    C = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            C[i,j] = cofactor(A, i, j)
    return C


def matrix_modinv(A: SquareMatrixInt, modulus: int) -> SquareMatrixModInt:
    C = cofactor_matrix(A) % modulus
    det_inv = integer_ring.modinv(det(A), modulus)

    return (det_inv * C.T) % modulus


def q_ary_lattice_basis(A: MatrixInt, modulus: int) -> SquareMatrixInt:
    m, n = A.shape
    assert n >= m
    # A = (A1 | A2)
    A1 = A[ : ,:m] # m x m
    A2 = A[ : ,m:] # m x (n - m)

    B11 = np.identity(m, dtype=int) # m x m
    B12 = (matrix_modinv(A1, modulus) @ A2) % modulus # m x (n - m)
    B21 = np.zeros((n - m, m), dtype=int)
    B22 = modulus * np.identity(n - m, dtype=int)

    print(B11)
    print()
    print(B12)
    print()
    print(B21)
    print()
    print(B22)

    return np.block([[B11, B12], [B21, B22]])