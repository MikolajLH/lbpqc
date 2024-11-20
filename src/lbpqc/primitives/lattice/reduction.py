import numpy as np


from lbpqc.primitives.lattice import full_rank_lattice as lattice


def gram_schmidt_process_coeff(basis: np.ndarray[float]) -> tuple[np.ndarray[float],np.ndarray[float]]:
    '''
    Performs Gram Schmidt orthogonalization of Vector Space basis.
    '''
    proj = lambda u,v: (np.dot(v,u) / np.dot(u,u))
    n,m = basis.shape
    B = basis.astype(float)

    new_basis = []
    coeff = np.zeros((n,m),float)
    for k in range(n):
        new_basis.append(B[k])
        coeff[k][0] = 1
        for j in range(k):
            coeff[k][j] = proj(new_basis[j], B[k])
            new_basis[k] -= coeff[k][j] * new_basis[j]
            
    return np.array(new_basis), coeff


def gram_schmidt_process(basis: np.ndarray[float]) -> np.ndarray[float]:
    return gram_schmidt_process_coeff(basis)[0]


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


def bounds(R: float, w: float, norm: float, alpha: float, all0: bool, eps: float = 1e-5) -> tuple[int, int]:
    K = np.sqrt(max(R**2 - w, 0)) / norm
    lower_bound = np.ceil(-K - alpha - eps)
    upper_bound = np.floor(K - alpha + eps)
    if all0:
        return 0, upper_bound
    return lower_bound, upper_bound


def enumerate(basis: np.ndarray, gs_basis: np.ndarray, coeff: np.ndarray, n:int, level: int,
              x: list[int], R: float, short_vec: np.ndarray, w: float, combination: list[int],
              eps: float = 1e-5) -> tuple[float, np.ndarray, list[int]]:
    if level < 0:
        new_vec = basis[0] * x[n - 1]
        for i in range(1, n):
            new_vec += basis[i] * x[n - 1 - i]
        new_R = np.sqrt(np.dot(new_vec, new_vec))
        if new_R > eps and R > new_R + eps:
            return new_R, new_vec, np.array(x)[::-1]
        return R, short_vec, combination
    else:
        new_R = R
        new_vec = np.copy(short_vec)
        new_combination = combination.copy()
        alpha = 0
        for i in range(level+1, n):
            alpha += x[n-1 - i] * coeff[i][level]
        all0 = True
        for i in range(len(x)):
            if x[i] != 0:
                all0 = False
                break
        lower_bound, upper_bound = bounds(new_R, w, np.sqrt(np.dot(gs_basis[level], gs_basis[level])), alpha, all0)
        i = upper_bound
        while i >= lower_bound:
            y = x.copy()
            y.append(i)
            res = enumerate(basis, gs_basis, coeff, n, level - 1, y, new_R, new_vec,
                            w + ((i + alpha)**2) * np.dot(gs_basis[level], gs_basis[level]), new_combination)
            if res[0] + eps < new_R:
                new_R = res[0]
                new_vec = res[1]
                new_combination = res[2]
                lower_bound, upper_bound = bounds(new_R, w, np.sqrt(np.dot(gs_basis[level], gs_basis[level])), alpha, all0)
            i -= 1
        return new_R, new_vec, new_combination


def SVP(basis: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[int]]:
    n, m = basis.shape
    gs_basis, coeff = gram_schmidt_process_coeff(basis)
    R = np.dot(basis[0], basis[0])
    short_vec = np.copy(basis[0])
    for i in range(1, n):
        new_R = np.dot(basis[i], basis[i])
        if new_R < R:
            R = new_R
            short_vec = np.copy(basis[i])
    new_R, new_vec, new_combination = enumerate(basis, gs_basis, coeff, n, n-1, [], R, short_vec, 0, [0]*n)
    return new_R, new_vec, basis, gs_basis, coeff, new_combination


def projected_lattice(basis: np.ndarray, start: int, end: int) -> np.ndarray:
    proj_basis = []
    gs_basis, coeff = gram_schmidt_process_coeff(basis)
    for i in range(start, end):
        tmp = np.copy(basis[i])
        for j in range(start):
            tmp -= gs_basis[j] * coeff[i][j]
        proj_basis.append(tmp)
    return np.array(proj_basis)


def lleaving_basis_vector(basis: np.ndarray, combination: list[int]) -> int:
    r = 0
    idx = -1
    for i in range(len(basis)):
        if combination[i] in [1, -1]:
            if r < np.dot(basis[i], basis[i]):
                idx = i
                r = np.dot(basis[i], basis[i])
    return idx


def BKZ(basis: np.ndarray, block_size: int) -> np.ndarray:
    n, m = basis.shape
    block_size = min(block_size, n)
    B = LLL(basis.astype(float))
    for i in range(n - block_size + 1):
        for j in range(block_size):
            proj_basis = projected_lattice(B, i + j, i + block_size)
            combination = SVP(proj_basis)[5]
            if np.dot(np.array(combination), np.array(combination)) > 1:
                new_basis = []
                for k in range(i + j):
                    new_basis.append(np.copy(B[k]))
                lifted_vec = [0]*m
                for k in range(len(combination)):
                    lifted_vec -= B[i+j+k] * combination[k]
                new_basis.append(lifted_vec)
                idx = lleaving_basis_vector(B, [0]*(i+j) + combination + [0]*(n-i-j-len(combination)))
                for k in range(i + j, n):
                    if k != idx:
                        new_basis.append(np.copy(B[k]))
                B = LLL(new_basis)
    return B


def HKZ(basis: np.ndarray) -> np.ndarray:
    n, m = basis.shape
    B = LLL(basis.astype(float))
    for i in range(n):
        proj_basis = projected_lattice(B, i, n)
        combination = SVP(proj_basis)[5]
        if np.dot(np.array(combination), np.array(combination)) > 1:
            new_basis = []
            for k in range(i):
                new_basis.append(np.copy(B[k]))
            lifted_vec = [0]*m
            for k in range(len(combination)):
                lifted_vec -= B[i+k] * combination[k]
            new_basis.append(lifted_vec)
            idx = lleaving_basis_vector(B, [0]*i + combination)
            for k in range(i, n):
                if k != idx:
                    new_basis.append(np.copy(B[k]))
            B = LLL(new_basis)
    return B