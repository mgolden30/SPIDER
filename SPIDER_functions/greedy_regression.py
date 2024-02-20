import numpy as np

def greedy_regression_pure_python(A, eps):
    """
    PURPOSE:
    The minimum of L = |A*c| in nested sparse subspaces, such that L
    increases minimally at each stage.

    INPUT:
    A - a matrix to look for sparse null vectors of

    OUTPUT:
    cs - columns of this matrix are the increasingly sparse approximate null
         vectors.
    residuals - vecnorm( A*cs );
    """
    m, n = A.shape

    fbounds = np.zeros((n, 2))
    cs = np.zeros((n, n))
    residuals = np.zeros(n)

    if m < n:
        print("error: matrix is underdetermined.")
        return None, None, None

    # first rotate A so it is square and upper triangular
    Q, R = np.linalg.qr(A)
    A = R[:n, :n]
    A0 = A.copy()

    I = np.ones(n, dtype=bool)  # logical vector indicating sparsity

    while n > 0:
        U, S, Vt = np.linalg.svd(A, full_matrices=False)
        V = Vt.T
        cs[I, n-1] = V[:, n-1]  # save out the smallest singular vector
        residuals[n-1] = S[n-1]

        if n == 1:
            break

        candidates = np.zeros(n)
        for i in range(n):
            a = A[:, i]
            alpha = 1 / np.linalg.norm(a)
            w = alpha * U.T @ a

            s = np.diag(S)
            bounds = [s[-1], s[-2]]

            f0 = lambda sigma: 1 - 1/alpha**2 * np.sum(w**2 / (s**2 - sigma**2 - eps**2))
            reg = lambda sigma: (s[-1]**2 - sigma**2 - eps**2) * (s[-2]**2 - sigma**2 - eps**2) * alpha**2 / (s[-1]**2 - s[-2]**2)
            f = lambda sigma: f0(sigma) * reg(sigma)

            maxit = 128
            threshold = 1e-130
            for j in range(maxit):
                g = np.sum(bounds) / 2  # bisection guess
                fg = f(g)
                if abs(fg) < threshold:
                    break

                if fg > 0:
                    bounds[0] = g
                else:
                    bounds[1] = g
            candidates[i] = g

        i_min = np.argmin(candidates)
        j = np.where(I)[0]
        I[j[i_min]] = False
        A = A0[:, I]
        n -= 1

    return cs, residuals, fbounds