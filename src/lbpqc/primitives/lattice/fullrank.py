from lbpqc.type_aliases import *


@enforce_type_check
def volume(lattice_basis: SquareMatrix) -> float:
    return abs(np.linalg.det(lattice_basis))


@enforce_type_check
def rank(lattice_basis: SquareMatrix) -> int:
    return lattice_basis.shape[0]


@enforce_type_check
def hadamard_ratio(lattice_basis: SquareMatrix) -> float:
    return (volume(lattice_basis) / np.linalg.norm(lattice_basis, axis=1).prod()) ** (1/rank(lattice_basis))


@enforce_type_check
def gaussian_expected_shortest_length(lattice_basis: SquareMatrix) -> float:
    n = rank(lattice_basis)
    return np.sqrt(n / (2 * np.pi * np.e)) * (volume(lattice_basis) ** (1/n))


@enforce_type_check
def transition_matrix(from_basis: SquareMatrix, to_basis: SquareMatrix) -> SquareMatrixInt:
    return np.rint(to_basis @ np.linalg.inv(from_basis)).astype(int)


@enforce_type_check
def babai_cvp(arbitrary_vector: Vector, lattice_basis: SquareMatrix) -> VectorInt:
    return np.rint(arbitrary_vector @ np.linalg.inv(lattice_basis)).astype(int)


@enforce_type_check
def dual_basis(lattice_basis: SquareMatrix) -> SquareMatrix:
    return np.linalg.inv(lattice_basis.T)
