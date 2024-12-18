from lbpqc.type_aliases import *
from lbpqc.primitives import matrix

def q_ary_basis(A: MatrixInt, q: int) -> MatrixModInt:
    m, n = A.shape
    B = np.block([[A.T], [q * np.identity(m, int)]])
    H, *_ = matrix.HNF(B)
    return H[:m, :m]


def dual_q_ary_basis(A: MatrixInt, q: int) -> MatrixModInt:
    m, n = A.shape
    Y = matrix.mod_left_kernel(A,q)
    return np.block([[Y], [ np.zeros((n, m - n), dtype=np.int64), q * np.identity(n, int)]])


def bai_galbraith_embedding(A: MatrixInt, b: VectorInt, q: int) -> SquareMatrixInt:
    #B = np.block([np.identity(m, int), A, -b.reshape(-1,1)])
    pass


def CVP_embedding(lattice_basis: SquareMatrix, v: Vector, M = 1) -> SquareMatrix:
    assert v.shape[0] == lattice_basis.shape[0]
    n = lattice_basis.shape[0]

    
    return np.block([[lattice_basis, np.zeros((n, 1), dtype=lattice_basis.dtype)], [v, np.array([M])]])


def subset_sum_lattice(sequence, S):
    n = len(sequence)
    M = np.identity(n + 1, dtype=float) * 2
    M[-1] = 1
    M[:-1, -1] = np.array(sequence, dtype=float)
    M[-1, -1] = S