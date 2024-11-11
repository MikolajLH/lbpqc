import numpy as np



def hnf(A: np.ndarray[int]):
    r'''
    Computes row-style Hermite Normal Form of a integer matrix A.
    '''
    assert A.dtype == int
    assert A.ndim == 2
    H = A.copy()
    m, n = H.shape
    p = min(m,n)
    k, j = 0, 0
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
        if H[k,j] < 0:
            H[k] = -H[k]
        # Reduce Rows
        b = H[k,j]
        for i in range(k+1, m):
            q = round(H[i,j] / b)
            H[i] -= q * H[k]
        # Check if column is done
        if np.all(H[k+1:, j] == 0):
            j += 1
            k += 1
            
    # Final reductions
    k = 0
    for j in range(p):
        if H[k,j] < 0:
            H[k] = -H[k]
        b = H[k,j]
        if b == 0: continue
        for i in range(k):
            q = H[i,j] // b
            H[i] -= q * H[k]
        k += 1
        
    return H