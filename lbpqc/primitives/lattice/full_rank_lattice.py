import numpy as np


def is_square_matrix(a: np.ndarray):
    r'''
    Predicate for checking if np.ndarray is square matrix,
    used mostly in various assertions across the codebase
    '''
    m, n, *_ = a.shape
    return len(_) == 0 and m == n


def lattice_det(lattice_basis: np.ndarray):
    r'''
    Full rank lattice determinant (covolume) is defined as an absolute value of determinant of a basis $B \subset \mathcal{R}^n $
    and it's an invariant independent of matrix B used to generate Lattice.

    $$ \det(L) = |\det(B)| $$
    '''
    assert is_square_matrix(lattice_basis)

    return abs(np.linalg.det(lattice_basis))


def lattice_dim(lattice_basis: np.ndarray):
    r'''
    Dimension of full rank lattice generated by basis $B \subset \mathcal{R}^n $  is equal to $n$.
    '''
    assert is_square_matrix(lattice_basis)

    n, _ = lattice_basis.shape
    return n


def hadamard_ratio(lattice_basis: np.ndarray) -> float:
    r'''
    Hadamard ratio is the indicator of the orthogonality of the lattice basis.

    The closer it is to 1 the more orthogonal lattice basis is.

    $$ \mathcal{H}(B) = (\frac{\det{L(B)}}{||b_1|| ||b_2|| \dotsb ||b_n|| })^{1/n} $$
    '''
    return (lattice_det(lattice_basis) / np.linalg.norm(lattice_basis, axis=1).prod()) ** (1/lattice_dim(lattice_basis))


def gaussian_expected_shortest_length(lattice_basis: np.ndarray) -> float:
    r'''
    computes the expected length of the shortest nonzero vector in a lattice using Gaussian formula
    '''
    n = lattice_dim(lattice_basis)
    return np.sqrt(n / (2 * np.pi * np.e)) * (lattice_det(lattice_basis) ** (1/n))


def lattice_change_of_basis_matrix(from_basis: np.ndarray, to_basis: np.ndarray) -> np.ndarray[int]:
    r'''
    Given two Lattice basis $ \mathcal{B} $ and $ \mathcal{B}^{*} $ computes unimodular matrix $ U $
    that represents transformation from one basis to another.
    '''
    assert is_square_matrix(from_basis) and is_square_matrix(to_basis)

    inv = np.linalg.inv(from_basis)
    return np.rint(np.array([w @ inv for w in to_basis])).astype(int)


def babai_closest_vector_algorithm(arbitrary_vector: np.ndarray, lattice_basis: np.ndarray) -> np.ndarray[int]:
    r'''
    Babai's algorithm to solve apprCVP.

    The more orthogonal lattice basis is, the better solution it provides.

    Let $ w \in \mathcal{R}^n $ be the vector for which we want to find the closest lattice vector $ v $ over the lattice basis $ B $.

    The algorithm works as follows:

    1. Write the vector $ w $  in coordinates of basis $ B $ by solving $ w = x \cdot B $.

    2. $ x \in \mathcal{R}^n $ so round it's coordinates to the nearest integers.

    3. return $ x $.
    '''
    return np.rint(arbitrary_vector @ np.linalg.inv(lattice_basis)).astype(int)


def knapsack_lattice_basis(int_sequence, S: int) -> np.ndarray:
    r'''
    Given knapsack problem defined by sequence $ \boldsymbol{r} = (r_1, r_2, ... , r_n) $ and sum $ S $
    computes associeted lattice basis

    TODO
    
    add latex formula for the matrix
    '''
    n = len(int_sequence)
    M = np.identity(n + 1, dtype=float) * 2
    M[-1] = 1
    M[:-1, -1] = np.array(int_sequence, dtype=float)
    M[-1, -1] = S
    return M



def NTRU_lattice_basis(NTRU_params, NTRU_public_key: np.ndarray[int]):
    r'''
    Given NTRU cryptosystem with parameters $ (N, p, q, d) $ and public key $ h $
    computes associeted lattice basis

    TODO

    add latex formula for the matrix
    '''
    N, p, q, d = NTRU_params
    h = NTRU_public_key

    M11 = np.identity(N, dtype=float)
    M21 = np.zeros((N,N),dtype=float)
    M22 = np.identity(N, dtype=float) * q
    M12 = np.array([np.roll(h, i) for i in range(N)], dtype=float)

    M_NTRU = np.block([[M11, M12], [M21, M22]])

    return M_NTRU