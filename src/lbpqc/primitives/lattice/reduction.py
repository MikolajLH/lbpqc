import numpy as np
from lbpqc.primitives.lattice import full_rank_lattice as lattice


def gram_schmidt_process(basis: np.ndarray[float]) -> np.ndarray[float]:
    '''
    Performs Gram Schmidt orthogonalization of Vector Space basis.
    '''
    proj = lambda u,v: (np.dot(v,u) / np.dot(u,u)) * u
    n,_ = basis.shape
    B = basis.astype(float)

    new_basis = []
    for k in range(n):
        new_basis.append(B[k])
        for j in range(k):
            new_basis[k] -= proj(new_basis[j], B[k])
            
    return np.array(new_basis)


def GLR_2dim(w1: np.ndarray, w2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    '''
    Gaussian Lattice reduction in dimension 2
    '''
    assert w1.shape == (2,)
    assert w2.shape == (2,)
    

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
    return v1, v2
    

def LLLdep(lattice_basis: np.ndarray, delta: float = 0.75):
    '''
    Implementation of LLL latice reduction algorithm
    '''
    n = lattice.lattice_dim(lattice_basis)
    B = lattice_basis.copy()
    Bstar = gram_schmidt_process(B)
    proj = lambda u,v: (np.dot(v,u) / np.dot(u,u))
    mu = lambda i, j: proj(Bstar[j], B[i])
    k = 1
    while k < n:
        for j in range(k-1, -1, -1):
            if abs(mu(k,j)) > 0.5:
                B[k] = B[k] - round(mu(k,j)) * B[j]
                Bstar =  gram_schmidt_process(B)
        if np.dot(Bstar[k], Bstar[k]) > (delta - (mu(k,k-1)**2)) * np.dot(Bstar[k-1], Bstar[k-1]):
           k = k + 1
        else:
            B[[k, k - 1]] = B[[k - 1, k]]
            Bstar =  gram_schmidt_process(B)
            k = max(k-1, 2)
    return B


# M = np.array([[19, 2, 32, 46, 3, 33], [15, 42, 11, 0, 3, 24], [43, 15, 0, 24, 4, 16], [20, 44, 44, 0, 18, 15], [0, 48, 35, 16, 31, 31], [48, 33, 32, 9, 1, 29]])
# B = np.array([[1, -1, 3], [1, 0, 5], [1, 2, 6]])

def LLL(lattice_basis: np.ndarray, delta: float = 0.75):
    n = lattice.lattice_dim(lattice_basis)
    B = lattice_basis.astype(float)
    while True:
        Bstar = gram_schmidt_process(B)
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


def babai_nearest_plane(lattice_basis: np.ndarray, w: np.ndarray):
    n = lattice.lattice_dim(lattice_basis)
    B = LLL(lattice_basis, 0.75)
    b = w.astype(float)
    for j in range(n - 1, -1, -1):
        Bstar = gram_schmidt_process(B)
        cj = round(np.dot(b, Bstar[j]) / np.dot(Bstar[j], Bstar[j]))
        b = b - cj * B[j]
    return w - b