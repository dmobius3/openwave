"""Saturated integral basis of the LEFT kernel {v in Z^n : v A = 0}.

Row-reduce A with a unimodular transform U so that U A = H is in row echelon
form. The rows of U whose H-row vanished are a basis of the left kernel, and it
is automatically SATURATED because U is unimodular, so those rows extend to a
Z-basis of Z^n. A rational RREF basis does not have this property, which is how
the rejected d3 came to generate a finite-index sublattice.
"""
def left_kernel(A):
    m, n = len(A), len(A[0])
    H = [row[:] for row in A]
    U = [[1 if i == j else 0 for j in range(m)] for i in range(m)]
    row = 0
    for col in range(n):
        # gcd-eliminate this column below `row` using integer row operations
        while True:
            nz = [r for r in range(row, m) if H[r][col]]
            if len(nz) <= 1: break
            nz.sort(key=lambda r: abs(H[r][col]))
            p = nz[0]
            for r in nz[1:]:
                q = H[r][col] // H[p][col]
                if q:
                    H[r] = [a - q*b for a, b in zip(H[r], H[p])]
                    U[r] = [a - q*b for a, b in zip(U[r], U[p])]
        nz = [r for r in range(row, m) if H[r][col]]
        if not nz: continue
        p = nz[0]
        H[row], H[p] = H[p], H[row]; U[row], U[p] = U[p], U[row]
        row += 1
        if row == m: break
    ker = [U[r] for r in range(m) if not any(H[r])]
    return ker, row
