from lbpqc.type_aliases import *


@enforce_type_check
def GSO(B: Matrix) -> Tuple[MatrixFloat, SquareMatrixFloat]:
    r'''

    Args:

    Returns:
    
    '''
    m, n = B.shape
    proj_coeff = lambda q, b: np.dot(b, q) / np.dot(q, q)
    B_star = B.astype(float)
    U = np.identity(m)

    for j in range(1, m):
        b = B_star[j].copy()
        for i in range(j):
            U[i,j] = proj_coeff(B_star[i], b)
            B_star[j] -= U[i][j] * B_star[i]
    
    # B = U.T @ B_star
    return B_star, U


def is_size_reduced(lattice_basis: Matrix) -> bool:
    r'''

    Args:

    Returns:
    
    '''
    _, U = GSO(lattice_basis)
    return np.all(np.abs(U[np.fromfunction(lambda i, j: i < j, U.shape).nonzero()]) <= 0.5)


def is_basis_vector_size_reduced(lattice_basis: Matrix, k: int) -> bool:
    _, U = GSO(lattice_basis)
    return np.all(np.abs(U[:k,k]) <= 0.5)


def lovasz_condition(lattice_basis: Matrix, delta: float) -> bool:
    r'''

    Args:

    Returns:
    
    '''
    norm2 = lambda x: np.sum(x * x, axis=1)
    G, U = GSO(lattice_basis)
    lhs = delta * norm2(G[:-1])
    rhs = norm2(G[1:] + np.diag(U, 1)[:, np.newaxis] * G[:-1])
    return np.all(lhs <= rhs)


def is_LLL_reduced(lattice_basis: Matrix, delta: float):
    r'''

    Args:

    Returns:
    
    '''
    return is_size_reduced(lattice_basis) and lovasz_condition(lattice_basis, delta)


def size_reduction_of_basis_vector(lattice_basis: Matrix, k: int):
    B = lattice_basis.astype(float)
    m, n = B.shape
    _, U = GSO(B)
    for j in range(k - 1, -1, -1):
        if abs(U[j, k]) > 0.5:
            B[k] -= np.rint(U[j,k]) * B[j]
            for i in range(m):
                U[i, k] -= round(U[j, k]) * U[i, j]
    return B, U

def size_reduction(lattice_basis: Matrix):
    B = lattice_basis.astype(float)
    m, n = B.shape
    _, U = GSO(B)

    for k in range(m - 1, -1, -1):
        for j in range(k - 1, -1, -1):
            B[k] -= round(U[j,k]) * B[j]
            for i in range(m):
                U[k, i] -= round(U[j, k]) * U[i, j]
    return B




def LLL(lattice_basis: SquareMatrix, delta: float = 0.75) -> SquareMatrixFloat:
    r'''

    Args:

    Returns:
    
    '''
    n = lattice_basis.shape[0]
    B = lattice_basis.astype(float)
    while True:
        Bstar, _ = GSO(B)
        # Reduction Step
        for i in range(1, n):
            for j in range(i-1, -1, -1):
                cij = round(np.dot(B[i], Bstar[j]) / np.dot(Bstar[j], Bstar[j]))
                B[i] = B[i] - cij * B[j]
        # Swap step
        exists = False
        for i in range(n - 1):
            u = np.dot(B[i + 1], Bstar[i]) / np.dot(Bstar[i], Bstar[i])
            r = u * Bstar[i] + Bstar[i + 1]
            if delta * np.dot(Bstar[i], Bstar[i]) > np.dot(r, r):
                B[[i, i + 1]] = B[[i + 1, i]]
                exists = True
                break
        if not exists:
            break
    return B


def babai_nearest_plane(lattice_basis: SquareMatrix, w: VectorFloat):
    r'''

    Args:

    Returns:
    
    '''
    n = lattice_basis.shape[0]
    B = LLL(lattice_basis, 0.75)
    b = w.astype(float)
    for j in range(n - 1, -1, -1):
        Bstar, _ = GSO(B)
        cj = round(np.dot(b, Bstar[j]) / np.dot(Bstar[j], Bstar[j]))
        b = b - cj * B[j]
    return w - b



def GLR_2dim(lattice_basis: SquareMatrix) -> SquareMatrixFloat:
    '''
    Gaussian Lattice reduction in dimension 2

    Args:

    Returns:
    
    '''
    if lattice_basis.shape != (2,2):
        raise ValueError(f"Lattice has to have rank 2 for gaussian reduction")
    
    w1 = lattice_basis[0]
    w2 = lattice_basis[1]

    v1 = w1.astype(float)
    v2 = w2.astype(float)
    if np.linalg.norm(v1) > np.linalg.norm(v2):
        v1, v2 = v2, v1

    while np.linalg.norm(v2) > np.linalg.norm(v1):
        m = round(np.dot(v1, v2) / np.dot(v1, v1))
        if m == 0:
            return v1, v2
        v2 = v2 - m * v1
        if np.linalg.norm(v1) > np.linalg.norm(v2):
            v1, v2 = v2, v1

    return np.array([v1, v2])