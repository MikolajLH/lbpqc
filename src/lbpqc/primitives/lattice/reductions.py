from lbpqc.type_aliases import *


@enforce_type_check
def GSO(B: Matrix) -> Tuple[MatrixFloat, SquareMatrixFloat]:
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
    _, U = GSO(lattice_basis)
    return np.all(np.abs(U[np.fromfunction(lambda i, j: i < j, U.shape).nonzero()]) <= 0.5)


def is_basis_vector_size_reduced(lattice_basis: Matrix, k: int) -> bool:
    _, U = GSO(lattice_basis)
    return np.all(np.abs(U[:k,k]) <= 0.5)


def lovasz_condition(lattice_basis: Matrix, delta: float) -> bool:
    norm2 = lambda x: np.sum(x * x, axis=1)
    G, U = GSO(lattice_basis)
    lhs = delta * norm2(G[:-1])
    rhs = norm2(G[1:] + np.diag(U, 1)[:, np.newaxis] * G[:-1])
    return np.all(lhs <= rhs)


def is_LLL_reduced(lattice_basis: Matrix, delta: float):
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



def _bounds(R: float, w: float, norm: float, alpha: float, all0: bool, eps: float = 1e-5) -> Tuple[int, int]:
    K = np.sqrt(max(R**2 - w, 0)) / norm
    lower_bound = np.ceil(-K - alpha - eps)
    upper_bound = np.floor(K - alpha + eps)
    if all0:
        return 0, upper_bound
    return lower_bound, upper_bound



def _enumerate(basis: SquareMatrixFloat, gs_basis: SquareMatrixFloat, coeff: np.ndarray, n:int, level: int,
              x: list[int], R: float, short_vec: np.ndarray, w: float, combination: list[int],
              eps: float = 1e-5) -> Tuple[float, VectorFloat, list[int]]:
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
            alpha += x[n-1 - i] * coeff[level][i]
        all0 = True
        for i in range(len(x)):
            if x[i] != 0:
                all0 = False
                break
        lower_bound, upper_bound = _bounds(new_R, w, np.sqrt(np.dot(gs_basis[level], gs_basis[level])), alpha, all0)
        i = upper_bound
        while i >= lower_bound:
            y = x.copy()
            y.append(i)
            res = _enumerate(basis, gs_basis, coeff, n, level - 1, y, new_R, new_vec,
                            w + ((i + alpha)**2) * np.dot(gs_basis[level], gs_basis[level]), new_combination)
            if res[0] + eps < new_R:
                new_R = res[0]
                new_vec = res[1]
                new_combination = res[2]
                lower_bound, upper_bound = _bounds(new_R, w, np.sqrt(np.dot(gs_basis[level], gs_basis[level])), alpha, all0)
            i -= 1
        return new_R, new_vec, new_combination


def SVP(basis: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[int]]:
    n, m = basis.shape
    gs_basis, coeff = GSO(basis)
    R = np.dot(basis[0], basis[0])
    short_vec = np.copy(basis[0])
    for i in range(1, n):
        new_R = np.dot(basis[i], basis[i])
        if new_R < R:
            R = new_R
            short_vec = np.copy(basis[i])
    new_R, new_vec, new_combination = _enumerate(basis, gs_basis, coeff, n, n-1, [], R, short_vec, 0, [0]*n)
    return new_R, new_vec, basis, gs_basis, coeff, new_combination


def projected_lattice(basis: SquareMatrix, start: int, end: int) -> SquareMatrixFloat:
    proj_basis = []
    gs_basis, coeff = GSO(basis)
    for i in range(start, end):
        tmp = np.copy(basis[i])
        for j in range(start):
            tmp -= gs_basis[j] * coeff[j][i]
        proj_basis.append(tmp)
    return np.array(proj_basis)


def lleaving_basis_vector(basis: SquareMatrix, combination: list[int]) -> int:
    r = 0
    idx = -1
    for i in range(len(basis)):
        if combination[i] in [1, -1]:
            if r < np.dot(basis[i], basis[i]):
                idx = i
                r = np.dot(basis[i], basis[i])
    return idx


def BKZ(basis: SquareMatrix, block_size: int) -> SquareMatrixFloat:
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


def HKZ(basis: SquareMatrix) -> SquareMatrixFloat:
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